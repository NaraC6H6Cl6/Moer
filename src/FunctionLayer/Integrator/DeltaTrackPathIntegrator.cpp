#include "DeltaTrackPathIntegrator.h"
#include "CoreLayer/Math/Warp.h"
#include "FastMath.h"
#include "FunctionLayer/Material/NullMaterial.h"

DeltaTrackPathIntegrator::DeltaTrackPathIntegrator(std::shared_ptr<Camera> _camera,
                                                   std::unique_ptr<Film> _film,
                                                   std::unique_ptr<TileGenerator> _tileGenerator,
                                                   std::shared_ptr<Sampler> _sampler,
                                                   int _spp,
                                                   int _renderThreadNum)
    : VolPathIntegrator(_camera,
                        std::move(_film),
                        std::move(_tileGenerator),
                        _sampler, _spp,
                        _renderThreadNum) {
}

Spectrum DeltaTrackPathIntegrator::Li(const Ray &initialRay,
                                      std::shared_ptr<Scene> scene) {
    Spectrum L{0.0};
    Spectrum T{1.0};
    Spectrum pathPdf{1.0};
    Spectrum lightPathPdf{1.0};
    const double eps = 1e-4;
    int nBounces = 0;
    bool isSpecularBounce = false;

    Ray ray = initialRay;
    std::shared_ptr<Medium> medium = nullptr;
    auto itsOpt = scene->intersect(ray);

    PathIntegratorLocalRecord evalLightRecord = evalEmittance(scene, itsOpt, ray);
    L += T * evalLightRecord.f;

    if (itsOpt.has_value()) {
        medium = getTargetMedium(itsOpt.value(), -ray.direction);
    } else {
        return L;
    }

    while (true) {
        MediumSampleRecord mRec;
        if (medium && medium->sampleDistance(&mRec, ray, itsOpt.value(), sampler->sample2D())) {
            if (T.isBlack()) {
                return L;
            }

            //* delta tracking

            Spectrum majorant = mRec.sigmaA + mRec.sigmaS + mRec.sigmaN;

            // eval emittance
            Spectrum Le = medium->evalEmittance(mRec.scatterPoint);
            if (!Le.isBlack()) {
                Spectrum P = T * mRec.tr * mRec.sigmaA * Le;
                Spectrum emitPdf = pathPdf * majorant * mRec.tr;
                L += P / emitPdf.average();
            }

            // sample medium event
            MediumSampleEvent event = medium->sampleEvent(mRec.scatterPoint, sampler->sample2D());
            if (event == MediumSampleEvent::a) {
                // absorption
                return L;
            } else if (event == MediumSampleEvent::s) {
                // real scattering
                nBounces++;
                if (nBounces > nPathLengthLimit) {
                    return L;
                }

                T *= mRec.sigmaS * mRec.tr;
                pathPdf *= mRec.sigmaS * mRec.tr;

                // sample luminaire
                Intersection mediumScatteringIts = fulfillScatteringPoint(mRec.scatterPoint, ray.direction, medium);
                for (int i = 0; i < nDirectLightSamples; ++i) {
                    PathIntegratorLocalRecord sampleLightRecord = sampleDirectLighting(scene, mediumScatteringIts, ray);
                    PathIntegratorLocalRecord evalScatterRecord = evalScatter(mediumScatteringIts, ray, sampleLightRecord.wi);
                    if (!sampleLightRecord.f.isBlack()) {
                        double misw = sampleLightRecord.isDelta ? 1.0 : MISWeight(sampleLightRecord.pdf, evalScatterRecord.pdf);
                        L += T * misw * sampleLightRecord.f * evalScatterRecord.f / sampleLightRecord.pdf / nDirectLightSamples / pathPdf.average();
                    }
                }

                // sample phase
                PathIntegratorLocalRecord sampleScatterRecord = sampleScatter(mediumScatteringIts, ray);
                if (sampleScatterRecord.f.isBlack()) {
                    return L;
                }
                T *= sampleScatterRecord.f;
                lightPathPdf = pathPdf;
                pathPdf *= sampleScatterRecord.pdf;

                // spawn new ray
                ray = Ray{mediumScatteringIts.position + sampleScatterRecord.wi * eps, sampleScatterRecord.wi};
                isSpecularBounce = false;
                itsOpt = scene->intersect(ray);
                auto [sampleIts, tr] = intersectIgnoreSurface(scene, ray, medium);
                auto evalLightRecord = evalEmittance(scene, sampleIts, ray);
                if (!evalLightRecord.f.isBlack()) {
                    double misw = sampleScatterRecord.isDelta ? 1.0 : MISWeight(sampleScatterRecord.pdf, evalLightRecord.pdf);
                    L += T * tr * evalLightRecord.f * misw / pathPdf.average();
                }
                continue;

            } else {
                // null scattering
                T *= mRec.sigmaN * mRec.tr;
                pathPdf *= mRec.sigmaN * mRec.tr;
                lightPathPdf *= majorant * mRec.tr;
            }

        } else {
            if (medium) {
                Spectrum majorant = mRec.sigmaA + mRec.sigmaS + mRec.sigmaN;
                T *= mRec.sigmaN * mRec.tr;
                pathPdf *= mRec.sigmaN * mRec.tr;
                lightPathPdf *= majorant * mRec.tr;
            }

            if (!itsOpt.has_value()) {
                // Accumulate contributions from infinite light sources
                PathIntegratorLocalRecord evalLightRecord = evalEmittance(scene, itsOpt, ray);
                if (!evalLightRecord.f.isBlack()) {
                    if (isSpecularBounce) {
                        L += T * evalLightRecord.f / pathPdf.average();
                    } else {
                        lightPathPdf *= evalLightRecord.pdf;
                        L += T * evalLightRecord.f / (pathPdf + lightPathPdf).average();
                    }
                }
                return L;
            }

            auto its = itsOpt.value();
            its.medium = medium;

            //* Direct Illumination
            for (int i = 0; i < nDirectLightSamples; ++i) {
                PathIntegratorLocalRecord sampleLightRecord = sampleDirectLighting(scene, its, ray);
                PathIntegratorLocalRecord evalScatterRecord = evalScatter(its, ray, sampleLightRecord.wi);

                if (!sampleLightRecord.f.isBlack()) {
                    //* Multiple importance sampling
                    double misw = sampleLightRecord.isDelta ? 1.0 : MISWeight(sampleLightRecord.pdf, evalScatterRecord.pdf);
                    L += T * misw * sampleLightRecord.f * evalScatterRecord.f / sampleLightRecord.pdf / nDirectLightSamples / pathPdf.average();
                }
            }

            if (its.material->getBxDF(its)->isNull()) {
                medium = getTargetMedium(its, ray.direction);
                ray = Ray{its.position + eps * ray.direction, ray.direction};
                itsOpt = scene->intersect(ray);
                if (!itsOpt) {
                    return L;
                }
                continue;
            }

            nBounces++;
            if (nBounces > nPathLengthLimit) {
                return L;
            }

            //* ----- BSDF Sampling -----
            PathIntegratorLocalRecord sampleScatterRecord = sampleScatter(its, ray);
            isSpecularBounce = sampleScatterRecord.isDelta;
            if (sampleScatterRecord.f.isBlack()) {
                return L;
            }
            T *= sampleScatterRecord.f;
            pathPdf *= sampleScatterRecord.pdf;

            //* Test whether the sampling ray hit the emitter

            ray = Ray{its.position + sampleScatterRecord.wi * eps, sampleScatterRecord.wi};
            itsOpt = scene->intersect(ray);

            auto [sampleIts, tr] = intersectIgnoreSurface(scene, ray, medium);

            auto evalLightRecord = evalEmittance(scene, sampleIts, ray);
            if (!evalLightRecord.f.isBlack()) {
                //* The sampling ray hit the emitter
                //* Multiple importance sampling
                double misw = sampleScatterRecord.isDelta ? 1.0 : MISWeight(sampleScatterRecord.pdf, evalLightRecord.pdf);
                L += T * tr * evalLightRecord.f * misw / pathPdf.average();
            }
        }
    }
    return L;
}

Spectrum DeltaTrackPathIntegrator::evalTransmittance(std::shared_ptr<Scene> scene,
                                                     const Intersection &its,
                                                     Point3d pointOnLight) const {

    float tmax = (pointOnLight - its.position).length();
    Ray shadowRay{its.position, normalize(pointOnLight - its.position), 1e-4, tmax - 1e-4};
    std::shared_ptr<Medium> medium = its.medium;
    Spectrum T(1);
    while (true) {
        // ratio tracking

        auto itsOpt = scene->intersect(shadowRay);

        if (medium) {
            if (!itsOpt) {
                Intersection itsTemp;
                itsTemp.medium = medium;
                itsTemp.position = pointOnLight;
                itsTemp.t = (pointOnLight - shadowRay.origin).length();
                MediumSampleRecord mRec;
                while (medium->sampleDistance(&mRec, shadowRay, itsTemp, sampler->sample2D())) {
                    // T *= (pN + tr * (pA + pS))
                    T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
                    shadowRay.origin = mRec.scatterPoint;
                    shadowRay.timeMax -= mRec.marchLength;
                    itsTemp.t -= mRec.marchLength;
                }
                T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
                return T;
            }

            if (!itsOpt->material->getBxDF(*itsOpt)->isNull()) {
                return 0.0;
            }

            MediumSampleRecord mRec;
            while (medium->sampleDistance(&mRec, shadowRay, itsOpt.value(), sampler->sample2D())) {
                T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
                shadowRay.origin = mRec.scatterPoint;
                shadowRay.timeMax -= mRec.marchLength;
                itsOpt.value().t -= mRec.marchLength;
            }

            T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
            medium = getTargetMedium(*itsOpt, shadowRay.direction);
            shadowRay.origin = itsOpt->position;
            shadowRay.timeMax -= itsOpt->t;
        } else {
            if (!itsOpt) {
                return T;
            }

            if (!itsOpt->material->getBxDF(*itsOpt)->isNull()) {
                return 0.0;
            }
            medium = getTargetMedium(*itsOpt, shadowRay.direction);
            shadowRay.origin = itsOpt->position;
            shadowRay.timeMax -= itsOpt->t;
        }
    }
}
