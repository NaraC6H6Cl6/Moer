﻿#include "Homogeneous.h"
#include "FastMath.h"

/// @brief Sample a scattering distance inside the (homogeneous) medium using naive but unbiased sampling strategy (i.e., inverse exponential sampling).
/// @param mRec Serve as return value.
/// @param ray Current ray inside medium.
/// @param its Current intersection on medium surface.
/// @param sample Two random number.
/// @return True if sampled a scattering point inside the medium, and thus fulfill marchLength, scatteringPoint, pdf, sigmaA, sigmaS, tr of mRec.
///         False if sampled a point on object boundary (surface), and thus fulfill marchLength, pdf, tr of mRec.
bool HomogeneousMedium::sampleDistance(MediumSampleRecord *mRec,
                                       const Ray &ray,
                                       const Intersection &its,
                                       Point2d sample) const {

    auto [x, y] = sample;

    // * randomly pick a channel/frequency and sample a distance.
    int channelIndex = int(x * nSpectrumSamples);
    double dist = -fm::log(1 - y) / mSigmaT[channelIndex];

    if (dist < its.t) {
        // * sampled a scattering point inside the medium.
        mRec->marchLength = dist;
        mRec->scatterPoint = ray.at(dist);
        mRec->sigmaA = mSigmaA;
        mRec->sigmaS = mSigmaS;
        mRec->sigmaN = Spectrum{0};
        mRec->tr = evalTransmittance(ray.origin, mRec->scatterPoint);
        // calculate pdf, i.e., sum of $1\n * \sigma_t^i e^{-\sigma_t^i * t}$.
        for (int i = 0; i < nSpectrumSamples; i++) {
            mRec->pdf[i] += mSigmaT[i] * fm::exp(-mSigmaT[i] * dist);
        }
        return true;
    } else {
        // sampled a point on object boundary (surface).
        mRec->marchLength = its.t;
        mRec->scatterPoint = ray.at(its.t);
        mRec->sigmaA = mSigmaA;
        mRec->sigmaS = mSigmaS;
        mRec->sigmaN = Spectrum{0};
        mRec->tr = evalTransmittance(ray.origin, its.position);
        // calculate discrete probility (instead of continuous probability density), i.e., sum of $1\n * e^{-\sigma_t^i * t_max}$.
        for (int i = 0; i < nSpectrumSamples; i++) {
            mRec->pdf[i] += fm::exp(-mSigmaT[i] * its.t);
        }
        return false;
    }
    // * Incomprehensible strategy that discrete probabilities and continuous probabilities share the same responsibility.
}

/// @brief Eval the exponential transmittance inside (homogeneous) medium. i.e., $e^{-\sigma_t * d}$ where $d$ stands for distance between two points inside medium.
/// @param from Start point.
/// @param dest End point.
/// @return Exponential transmittance per channel or frequency.
Spectrum HomogeneousMedium::evalTransmittance(Point3d from,
                                              Point3d dest) const {
    double dist = (dest - from).length();
    Spectrum tr;
    for (int i = 0; i < nSpectrumSamples; i++) {
        tr[i] = fm::exp(-mSigmaT[i] * dist);
    }
    return tr;
}

MediumSampleEvent HomogeneousMedium::sampleEvent(Point3d pos,
                                                 Point2d sample) const {
    auto [x, y] = sample;
    int channelIndex = int(x * nSpectrumSamples);

    if (y < mSigmaA[channelIndex] / mSigmaT[channelIndex]) {
        return MediumSampleEvent::a;
    } else {
        return MediumSampleEvent::s;
    }
}

Spectrum HomogeneousMedium::evalEmittance(Point3d pos) const {
    // TODO
    return Spectrum{.0};
}
