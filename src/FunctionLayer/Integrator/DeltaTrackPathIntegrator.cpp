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
    Spectrum neePdf{1.0};
    const double eps = 1e-4;
    int nBounces = 0;
    bool isSpecularBounce = false;
    bool atMediumBoundary = false;

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

        if (T.isBlack()) {
            break;
        }

        MediumSampleRecord mRec;
        if (medium && !atMediumBoundary) {
            atMediumBoundary = !medium->sampleDistance(&mRec, ray, itsOpt.value(), sampler->sample2D());

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
                neePdf = pathPdf;
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
                nBounces++;
                if (nBounces > nPathLengthLimit || !itsOpt.has_value()) {
                    break;
                }
                continue;

            } else {
                // null scattering
                ray.origin = mRec.scatterPoint + ray.direction * eps;
                ray.timeMax -= mRec.marchLength;
                itsOpt.value().t -= mRec.marchLength;
                T *= mRec.sigmaN * mRec.tr;
                pathPdf *= mRec.sigmaN * mRec.tr;
                neePdf *= majorant * mRec.tr;
                continue;
            }

        } else {

            auto its = itsOpt.value();
            its.medium = medium;
            atMediumBoundary = false;

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
                    break;
                }
                continue;
            }

            //* ----- BSDF Sampling -----
            PathIntegratorLocalRecord sampleScatterRecord = sampleScatter(its, ray);
            isSpecularBounce = sampleScatterRecord.isDelta;
            if (sampleScatterRecord.f.isBlack()) {
                break;
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

            nBounces++;
            if (nBounces > nPathLengthLimit || !itsOpt.has_value()) {
                break;
            }
        }
    }
    // Accumulate contributions from infinite light sources
    evalLightRecord = evalEnvLights(scene, ray);
    if (!evalLightRecord.f.isBlack()) {
        if (isSpecularBounce) {
            L += T * evalLightRecord.f / pathPdf.average();
        } else {
            neePdf *= evalLightRecord.pdf;
            L += T * evalLightRecord.f / (pathPdf + neePdf).average();
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
    return T;
}

std::pair<std::optional<Intersection>, Spectrum>
DeltaTrackPathIntegrator::intersectIgnoreSurface(std::shared_ptr<Scene> scene,
                                                 const Ray &ray,
                                                 std::shared_ptr<Medium> medium) const {

    const double eps = 1e-5;
    Vec3d dir = ray.direction;

    Spectrum T{1.0};
    Ray marchRay{ray.origin + dir * eps, dir};
    std::shared_ptr<Medium> currentMedium = medium;

    Point3d lastScatteringPoint = ray.origin;
    auto testRayItsOpt = scene->intersect(marchRay);

    // calculate the transmittance of last segment from lastScatteringPoint to testRayItsOpt.
    while (true) {

        // corner case: infinite medium or infinite light source.
        if (!testRayItsOpt.has_value()) {
            if (currentMedium != nullptr) {
                T = Spectrum(0.0);
            }
            return {testRayItsOpt, T};
        }

        auto testRayIts = testRayItsOpt.value();

        // corner case: non-null surface
        if (testRayIts.material != nullptr) {
            if (!testRayIts.material->getBxDF(testRayIts)->isNull()) {
                if (currentMedium != nullptr) {
                    // ratio tracking
                    Intersection itsTemp = testRayIts;
                    MediumSampleRecord mRec;
                    while (medium->sampleDistance(&mRec, marchRay, itsTemp, sampler->sample2D())) {
                        // T *= (pN + tr * (pA + pS))
                        T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
                        marchRay.origin = mRec.scatterPoint;
                        marchRay.timeMax -= mRec.marchLength;
                        itsTemp.t -= mRec.marchLength;
                    }
                    T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
                }
                return {testRayItsOpt, T};
            }
        }

        // hit a null surface, calculate T
        if (currentMedium != nullptr) {
            // ratio tracking
            Intersection itsTemp = testRayIts;
            MediumSampleRecord mRec;
            while (medium->sampleDistance(&mRec, marchRay, itsTemp, sampler->sample2D())) {
                // T *= (pN + tr * (pA + pS))
                T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
                marchRay.origin = mRec.scatterPoint;
                marchRay.timeMax -= mRec.marchLength;
                itsTemp.t -= mRec.marchLength;
            }
            T *= (mRec.sigmaN + mRec.tr * (mRec.sigmaA + mRec.sigmaS)) / (mRec.sigmaA + mRec.sigmaS + mRec.sigmaN);
        }

        // update medium
        currentMedium = getTargetMedium(testRayIts, dir);

        // update ray and intersection point.
        marchRay.origin = testRayIts.position + dir * eps;
        lastScatteringPoint = testRayIts.position;
        testRayItsOpt = scene->intersect(marchRay);
    }
    return {testRayItsOpt, T};
}
