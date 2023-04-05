/**
 * @file VolPathIntegrator.cpp
 * @author Chenxi Zhou
 * @brief 
 * @version 0.1
 * @date 2022-09-22
 * 
 * @copyright NJUMeta (c) 2022 
 * www.njumeta.com
 */

#include "VolPathIntegrator.h"
#include "CoreLayer/Math/Warp.h"
#include "FastMath.h"
#include "FunctionLayer/Material/NullMaterial.h"

VolPathIntegrator::VolPathIntegrator(std::shared_ptr<Camera> _camera,
                                     std::unique_ptr<Film> _film,
                                     std::unique_ptr<TileGenerator> _tileGenerator,
                                     std::shared_ptr<Sampler> _sampler,
                                     int _spp,
                                     int _threadThreadNum)
    : AbstractPathIntegrator(_camera, std::move(_film),
                             std::move(_tileGenerator),
                             _sampler, _spp,
                             _threadThreadNum)
{

}

Spectrum VolPathIntegrator::Li(const Ray &initialRay, std::shared_ptr<Scene> scene) {
    const double eps = 1e-4;
    Spectrum L{.0};
    Spectrum throughput{1.0};
    Ray ray = initialRay;
    int nBounces = 0;

    auto itsOpt = scene->intersect(ray);
    std::shared_ptr<Medium> medium = nullptr;

    MediumSampleRecord mRec;
    while(true) {

        //* All rays generated by bsdf/phase sampling which might hold radiance
        //* will be counted at the end of the loop
        //* except the ray generated by camera which is considered here.
        //* The radiance of environment map will be counted while no intersection is found but not here.
        if (nBounces == 0) {
            if (itsOpt.has_value()) {
                PathIntegratorLocalRecord evalLightRecord = evalEmittance(scene, itsOpt, ray);
                L += throughput * evalLightRecord.f;
            }    
        }

        //* No intersection. Add possible radiance from environment map.
        if (!itsOpt.has_value()){
            auto envMapRecord=evalEmittance(scene,itsOpt,ray);
            L += throughput * envMapRecord.f;
            break;
        }
        
        auto its = itsOpt.value();

        nBounces++;

        // * Ignore null materials using dynamic_cast.
        if(dynamic_cast<NullMaterial*>(its.material.get())!=nullptr){
            nBounces--;
            // ray should be immersed in possible medium.
            ray = Ray{its.position + ray.direction * eps, ray.direction};
            itsOpt=scene->intersect(ray);
            continue;
        }

        nBounces++;
        double pSurvive = russianRoulette(scene, itsOpt.value(), throughput, nBounces);
        if (randFloat() > pSurvive)
            break;
        throughput /= pSurvive;

        // * Ray currently travel inside medium and will continue.
    //     if (medium &&
    //         medium->sampleDistance(&mRec, ray, itsOpt.value(), sampler->sample2D()))
    //     {
    //         throughput *= mRec.tr * mRec.sigmaS / mRec.pdf;
    //         //* ----- Luminaire Sampling -----
    //         for (int i = 0; i < nDirectLightSamples; ++i) {
    //             PathIntegratorLocalRecord sampleLightRecord = sampleDirectLighting(scene, mRec, ray, medium);
    //             PathIntegratorLocalRecord evalScatterRecord = evalScatter(scene, mRec, ray, sampleLightRecord.wi, medium);
            
    //             if (!sampleLightRecord.f.isBlack()) {
    //                 double misw = MISWeight(sampleLightRecord.pdf, evalScatterRecord.pdf);
    //                 if (sampleLightRecord.isDelta)
    //                     misw = 1.0;
    //                 L += throughput * sampleLightRecord.f * evalScatterRecord.f
    //                     / sampleLightRecord.pdf * misw
    //                     / nDirectLightSamples;
    //             }
    //         }
    //         //* ----- Phase Sampling -----
    //         PathIntegratorLocalRecord sampleScatterRecord = sampleScatter(scene, mRec, ray, medium);
    //         if (sampleScatterRecord.f.isBlack())
    //             break;
    //         throughput *= sampleScatterRecord.f / sampleScatterRecord.pdf;
    //         ray = Ray{mRec.scatterPoint + sampleScatterRecord.wi * 1e-4, sampleScatterRecord.wi};
    //         itsOpt = scene->intersect(ray);

    //         auto[sampleIts, tr]
    //             = intersectIgnoreSurface(scene, ray, medium);
    //         auto evalLightRecord = evalEmittance(scene, sampleIts, ray);
    //         if (!evalLightRecord.f.isBlack()) {
    //             double misw = MISWeight(sampleScatterRecord.pdf, evalLightRecord.pdf);
    //             if (sampleScatterRecord.isDelta)
    //                 misw = 1.0;

    //             L += throughput * tr * evalLightRecord.f * misw;
    //         }

    //         if (!itsOpt.has_value())
    //             break;

    //     } 
    //     else {
    //         // * Ray currently travel inside medium, but will flee from medium.
    //         if (medium) {
    //             throughput *= mRec.tr / mRec.pdf;
    //         }
            

    //         //* Evaluate the radiance only when the ray is generated by camera
    //         if (nBounces == 0) {
    //             if (!evalLightRecord.f.isBlack()) {
    //                 L += throughput * evalLightRecord.f;
    //             }
    //         }
    //         if (!itsOpt.has_value())
    //             break;

    //         auto its = itsOpt.value();
    //         //* ----- Handle special surface -----
    //         if (its.material->type & EMaterialType::Null) {
    //             medium = getTargetMedium(ray, its, ray.direction);
    //             ray = Ray{its.position + 1e-4 * ray.direction, ray.direction};
    //             itsOpt = scene->intersect(ray);
    //             evalLightRecord = evalEmittance(scene, itsOpt, ray);
    //             continue;
    //         }

    //         //* ----- Direct Illumination -----
    //         for (int i = 0; i < nDirectLightSamples; ++i) {
    //             PathIntegratorLocalRecord sampleLightRecord = sampleDirectLighting(scene, its, ray);
    //             PathIntegratorLocalRecord evalScatterRecord = evalScatter(scene, its, ray, sampleLightRecord.wi);

    //             if (!sampleLightRecord.f.isBlack()) {
    //                 //* Multiple importance sampling
    //                 double misw = MISWeight(sampleLightRecord.pdf, evalScatterRecord.pdf);
    //                 if (sampleLightRecord.isDelta)
    //                     misw = 1.0;
    //                 L += throughput * sampleLightRecord.f * evalScatterRecord.f
    //                      / sampleLightRecord.pdf * misw
    //                      / nDirectLightSamples;
    //             }
    //         }

    //         //* ----- BSDF Sampling -----
    //         PathIntegratorLocalRecord sampleScatterRecord = sampleScatter(scene, its, ray);
    //         if (sampleScatterRecord.f.isBlack())
    //             break;
    //         throughput *= sampleScatterRecord.f / sampleScatterRecord.pdf;

    //         medium = getTargetMedium(ray, its, sampleScatterRecord.wi);

    //         //* Test whether the sampling ray hit the emitter
    //         const double eps = 1e-4;
    //         ray = Ray{its.position + sampleScatterRecord.wi * eps, sampleScatterRecord.wi};
    //         itsOpt = scene->intersect(ray);
          
    //         auto [sampleIts, tr] 
    //             = intersectIgnoreSurface(scene, ray, medium);

    //         evalLightRecord = evalEmittance(scene, sampleIts, ray);
    //         if (!evalLightRecord.f.isBlack()) {
    //             //* The sampling ray hit the emitter
    //             //* Multiple importance sampling

    //             double misw = MISWeight(sampleScatterRecord.pdf, evalLightRecord.pdf);
    //             if (sampleScatterRecord.isDelta)
    //                 misw = 1.0;

    //             L += throughput * tr * evalLightRecord.f * misw;
    //         }

    //         //* Terminate if ray escape the scene
    //         if (!itsOpt.has_value())
    //             break;

    //     }
    }

    return L;
}

/// @brief Eval surface or infinite light source emittance and take medium transmittance into account.
/// @param scene Scene description. Multiple shadow ray intersections will be performed.
/// @param itsOpt Current intersection point. If there's no intersection, eval the radiance of environment light.
/// @param ray Current ray.
/// @return current ray direction, obtained light radiance and corresponding solid angle dependent pdf.
PathIntegratorLocalRecord VolPathIntegrator::evalEmittance(std::shared_ptr<Scene> scene, 
                                                           std::optional<Intersection> itsOpt, 
                                                           const Ray &ray)
{
    Vec3d wo = -ray.direction;
    Spectrum LEmission(0.0);
    double pdfDirect = 1.0;
    if (!itsOpt.has_value())
    {
        auto record = evalEnvLights(scene, ray);
        LEmission = record.f;
        pdfDirect = record.pdf;
    }
    else if (itsOpt.value().object && itsOpt.value().object->getLight())
    {
        auto its = itsOpt.value();
        Normal3d n = its.geometryNormal;
        auto light = itsOpt.value().object->getLight();
        auto record = light->eval(ray, its, ray.direction);
        LEmission = record.s;
        Intersection tmpIts;
        tmpIts.position = ray.origin;
        pdfDirect = record.pdfDirect * chooseOneLightPdf(scene, tmpIts, ray, light);
    }
    Spectrum transmittance(1.0);
    return {ray.direction, transmittance * LEmission, pdfDirect, false}; 

}

PathIntegratorLocalRecord VolPathIntegrator::sampleDirectLighting(std::shared_ptr<Scene> scene, 
                                                                  const Intersection &its, 
                                                                  const Ray &ray)
{
    auto [light, pdfChooseLight] = chooseOneLight(scene, its, ray, sampler->sample1D());
    auto record = light->sampleDirect(its, sampler->sample2D(), ray.timeMin);
    double pdfDirect = record.pdfDirect * pdfChooseLight; // pdfScatter with respect to solid angle
    Vec3d dirScatter = record.wi;
    Spectrum Li = record.s;
    Point3d posL = record.dst;
    Point3d posS = its.position;
    Ray visibilityTestingRay(posL - dirScatter * 1e-4, -dirScatter, ray.timeMin, ray.timeMax);
    auto [visibilityTestingIts, tr] 
        = intersectIgnoreSurface(scene, visibilityTestingRay, nullptr);
    if (!visibilityTestingIts.has_value() || visibilityTestingIts->object != its.object || (visibilityTestingIts->position - posS).length2() > 1e-6)
    {
        tr = 0.0;
    }
    return {dirScatter, Li * tr, pdfDirect, record.isDeltaPos};

}

PathIntegratorLocalRecord VolPathIntegrator::evalScatter(std::shared_ptr<Scene> scene,
                                                      const Intersection &its,
                                                      const Ray &ray,
                                                      const Vec3d &dirScatter)
{
    if (its.material != nullptr)
    {
        std::shared_ptr<BxDF> bxdf = its.material->getBxDF(its);
        Normal3d n = its.geometryNormal;
        double wiDotN = fm::abs(dot(n, dirScatter));
        Vec3d wi = its.toLocal(dirScatter);
        Vec3d wo = its.toLocal(-ray.direction);
        return {
            dirScatter,
            bxdf->f(wo, wi,false) * wiDotN,
            bxdf->pdf(wo, wi),
            false};
    }
    else
    {
        // todo: eval phase function
        return {};
    }
}

PathIntegratorLocalRecord VolPathIntegrator::sampleScatter(std::shared_ptr<Scene> scene,
                                                        const Intersection &its,
                                                        const Ray &ray)
{
    if (its.material != nullptr)
    {
        Vec3d wo = its.toLocal(-ray.direction);
        std::shared_ptr<BxDF> bxdf = its.material->getBxDF(its);
        Vec3d n = its.geometryNormal;
        BxDFSampleResult bsdfSample = bxdf->sample(wo, sampler->sample2D(),false);
        double pdf = bsdfSample.pdf;
        Vec3d dirScatter = its.toWorld(bsdfSample.directionIn);
        double wiDotN = fm::abs(dot(dirScatter, n));
        return {dirScatter, bsdfSample.s * wiDotN, pdf, BxDF::MatchFlags(bsdfSample.bxdfSampleType,BXDF_SPECULAR)};
    }
    else
    {
        // todo: sample phase function
        return {};
    }
}

double VolPathIntegrator::russianRoulette(std::shared_ptr<Scene> scene,
                                       const Intersection &its,
                                       const Spectrum &throughput,
                                       int nBounce)
{
    // double pSurvive = std::min(pRussianRoulette, throughput.sum());
    double pSurvive = pRussianRoulette;
    if (nBounce > nPathLengthLimit)
        pSurvive = 0.0;
    if (nBounce <= 2)
        pSurvive = 1.0;
    return pSurvive;
}

std::pair<std::shared_ptr<Light>, double> 
VolPathIntegrator::chooseOneLight(std::shared_ptr<Scene> scene,
                                  const Intersection &its,
                                  const Ray &ray,
                                  double lightSample)
{
    // uniformly weighted
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    int numLights = lights->size();
    int lightID = std::min(numLights - 1, (int)(lightSample * numLights));
    std::shared_ptr<Light> light = lights->operator[](lightID);
    return {light, 1.0 / numLights};
}

double VolPathIntegrator::chooseOneLightPdf(std::shared_ptr<Scene> scene,
                                            const Intersection &its,
                                            const Ray &ray,
                                            std::shared_ptr<Light> light)
{
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    int numLights = lights->size();
    return 1.0 / numLights;
}

PathIntegratorLocalRecord VolPathIntegrator::evalEnvLights(std::shared_ptr<Scene> scene,
                                                           const Ray &ray)
{
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    Spectrum L(0.0);
    double pdf = 0.0;
    for (auto light : *lights)
    {
        auto record = light->evalEnvironment(ray);
        L += record.s;
        pdf += record.pdfEmitDir;
    }
    return {-ray.direction, L, pdf};
}

std::shared_ptr<Medium> VolPathIntegrator::getTargetMedium(const Ray &ray, 
                                                           const Intersection &its,
                                                           Vec3d wi) const 
{
    bool scatterToOutSide = dot(its.geometryNormal, wi) > 0;
    bool scatterToinSide = dot(its.geometryNormal, wi) < 0;
    if (scatterToOutSide) 
        return its.material->getOutsideMedium();
    else if (scatterToinSide)
        return its.material->getInsideMedium();

    return nullptr;
}

std::pair<std::optional<Intersection>, Spectrum> 
VolPathIntegrator::intersectIgnoreSurface(std::shared_ptr<Scene> scene, 
                                          const Ray &ray,
                                          std::shared_ptr<Medium> medium) const
{
    //TODO add an argument indicates which surface to ignore
    //! Infinity medium is not considered
    Ray marchingRay = ray;
    std::shared_ptr<Medium> currentMedium = medium;
    
    auto its = scene->intersect(marchingRay);
    Spectrum tr(1.0);

    while(its.has_value()) {
        tr *= (currentMedium ? currentMedium->evalTransmittance(ray.origin, its->position) : 1);

        if (its->material->type & EMaterialType::Null) {
            //* Continue the ray when the surface is ignored
            const double eps = 1e-4;
            marchingRay = Ray {its->position + eps * marchingRay.direction, marchingRay.direction};
            currentMedium = getTargetMedium(marchingRay, its.value(), marchingRay.direction);
        } else {
            //* Intersect on surface which can't be ignored
            break;
        }
        its = scene->intersect(marchingRay);
    }
    return {its, tr};
}

PathIntegratorLocalRecord 
VolPathIntegrator::sampleDirectLighting(std::shared_ptr<Scene> scene,
                                        const MediumSampleRecord &mRec,
                                        const Ray &ray,
                                        std::shared_ptr<Medium> medium) const 
{
    auto [light, pdfChooseLight]
        = chooseOneLight(scene, mRec, ray, sampler->sample1D());

    auto record 
        = light->sampleDirect(mRec, sampler->sample2D(), ray.timeMin);

    double pdfDirect = record.pdfDirect * pdfChooseLight;
    Vec3d dirScatter = record.wi;
    Spectrum Li = record.s;
    Point3d posL = record.dst,
            PosS = mRec.scatterPoint;

    Ray shadowRay{mRec.scatterPoint + dirScatter * 1e-4, dirScatter};
    auto [shadowIts, tr] 
        = intersectIgnoreSurface(scene, shadowRay, medium);
    if (!shadowIts.has_value() || (shadowIts->position - record.dst).length2() > 1e-6 ) {
        tr = 0.0;
    }

    return {dirScatter, Li * tr, pdfDirect, record.isDeltaPos};

}


std::pair<std::shared_ptr<Light>, double>
VolPathIntegrator::chooseOneLight(std::shared_ptr<Scene> scene,
                                  const MediumSampleRecord &mRec,
                                  const Ray &ray,
                                  double lightSample) const 
{
    // uniformly weighted
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    int numLights = lights->size();
    int lightID = std::min(numLights - 1, (int)(lightSample * numLights));
    std::shared_ptr<Light> light = lights->operator[](lightID);
    return {light, 1.0 / numLights};

}

PathIntegratorLocalRecord
VolPathIntegrator::evalScatter(std::shared_ptr<Scene> scene,
                               const MediumSampleRecord &mRec,
                               const Ray &ray,
                               const Vec3d &wi,
                               std::shared_ptr<Medium> medium) const
{
    auto [phaseValue, phasePdf, isDelta]
        = medium->evalPhase(-ray.direction, wi, mRec.scatterPoint);
    return {wi, phaseValue, phasePdf, isDelta};
}

PathIntegratorLocalRecord 
VolPathIntegrator::sampleScatter(std::shared_ptr<Scene> scene,
                                 const MediumSampleRecord &mRec,
                                 const Ray &ray,
                                 std::shared_ptr<Medium> medium) const
{
    auto [wi, phaseValue, phasePdf, isDelta]
        = medium->samplePhase(-ray.direction, mRec.scatterPoint, sampler->sample2D());
    return {wi, phaseValue, phasePdf, isDelta};
}  
