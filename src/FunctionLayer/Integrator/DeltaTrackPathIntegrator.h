#pragma once

#include "VolPathIntegrator.h"
#include "CoreLayer/Ray/Ray.h"
#include "AbstractPathIntegrator.h"
#include "FunctionLayer/Medium/Medium.h"

/**
 * @brief Bidirectional path-tracing integrator which
 * take volume into consideration
 * @ingroup Integrator
 */

class DeltaTrackPathIntegrator : public VolPathIntegrator {

public:
    DeltaTrackPathIntegrator(std::shared_ptr<Camera> _camera,
                             std::unique_ptr<Film> _film,
                             std::unique_ptr<TileGenerator> _tileGenerator,
                             std::shared_ptr<Sampler> _sampler,
                             int _spp,
                             int _renderThreadNum = 4);

    virtual Spectrum Li(const Ray &ray,
                        std::shared_ptr<Scene> scene) override;

    // virtual PathIntegratorLocalRecord evalEmittance(std::shared_ptr<Scene> scene,
    //                                                std::optional<Intersection> itsOpt,
    //                                                const Ray &ray) override;

    // virtual PathIntegratorLocalRecord sampleDirectLighting(std::shared_ptr<Scene> scene,
    //                                                        const Intersection &its,
    //                                                        const Ray &ray) override;

    // virtual PathIntegratorLocalRecord evalScatter(const Intersection &its,
    //                                               const Ray &ray,
    //                                               const Vec3d &wi) override;

    // virtual PathIntegratorLocalRecord sampleScatter(const Intersection &its,
    //                                                 const Ray &ray) override;

    // virtual double russianRoulette(const Spectrum &T,
    //                                int nBounce) override;

    // virtual std::pair<std::shared_ptr<Light>, double> chooseOneLight(std::shared_ptr<Scene> scene,
    //                                                                  double lightSample);

    // virtual double chooseOneLightPdf(std::shared_ptr<Scene> scene,
    //                                  std::shared_ptr<Light> light);

    // virtual PathIntegratorLocalRecord evalEnvLights(std::shared_ptr<Scene> scene,
    //                                                 const Ray &ray);

    virtual Spectrum evalTransmittance(std::shared_ptr<Scene> scene,
                                       const Intersection &its,
                                       Point3d pointOnLight) const;
};
