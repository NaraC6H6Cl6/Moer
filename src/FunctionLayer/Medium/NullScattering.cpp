#include "NullScattering.h"
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>
#include <CoreLayer/Adapter/Random.h>
#include <FunctionLayer/Sampler/Independent.h>
#include <ResourceLayer/File/FileUtils.h>

using FloatGridSampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 0, false>;
using Vec3fGridSampler = nanovdb::SampleFromVoxels<nanovdb::Vec3fGrid::TreeType, 0, false>;

NullScatteringMedium::NullScatteringMedium(std::string _gridFilePath,
                                           std::shared_ptr<PhaseFunction> _phase,
                                           TransformMatrix3D _transformMatrix,
                                           Spectrum _sigma_a,
                                           Spectrum _sigma_s)
    : Medium(_phase),
      transformMatrix(_transformMatrix),
      invTransformMatrix(_transformMatrix.getInverse()),
      sigma_a(_sigma_a),
      sigma_s(_sigma_s) {

    std::string fullGridFilePath = FileUtils::getWorkingDir() + _gridFilePath;
    densityGridBuffer = nanovdb::io::readGrid(fullGridFilePath, "density", 1);
    if (!densityGridBuffer) {
        std::cout << ".nvdb file must contains density grid\n";
        exit(1);
    }
    densityGrid = densityGridBuffer.grid<float>();

    emissionGridBuffer = nanovdb::io::readGrid(fullGridFilePath, "emission", 1);
    emissionGrid = emissionGridBuffer.grid<nanovdb::Vec3f>();

    // Assume every voxel is cube
    voxelSize = densityGrid->voxelSize()[0];

    auto bbox = densityGrid->worldBBox();
    Point3d boxMin = Point3d(bbox.min()[0], bbox.min()[1], bbox.min()[2]);
    Point3d boxMax = Point3d(bbox.max()[0], bbox.max()[1], bbox.max()[2]);

    //* Compute density grid bound
    minIndex[0] = densityGrid->indexBBox().min().x();
    minIndex[1] = densityGrid->indexBBox().min().y();
    minIndex[2] = densityGrid->indexBBox().min().z();

    maxIndex[0] = densityGrid->indexBBox().max().x();
    maxIndex[1] = densityGrid->indexBBox().max().y();
    maxIndex[2] = densityGrid->indexBBox().max().z();

    Spectrum sigma_t = sigma_a + sigma_s;
    float minVal, maxVal;
    densityGrid->tree().extrema(minVal, maxVal);
    majorant[0] = maxVal * sigma_t[0];
    majorant[1] = maxVal * sigma_t[1];
    majorant[2] = maxVal * sigma_t[2];
}

bool NullScatteringMedium::sampleDistance(MediumSampleRecord *mRec,
                                          const Ray &ray,
                                          const Intersection &its,
                                          Point2d sample) const {

    auto [x, y] = sample;

    // * randomly pick a channel/frequency and sample a distance.
    int channelIndex = int(x * nSpectrumSamples);
    double distance = -fm::log(1 - y) / (majorant[channelIndex]);

    if (distance < its.t) {
        // * sampled a scattering point inside the medium.
        mRec->marchLength = distance;
        mRec->scatterPoint = ray.at(distance);

        Point3d index = worldToIndex(mRec->scatterPoint);
        float density = sampleFromGrid(index, densityGrid);

        mRec->sigmaA = density * sigma_a;
        mRec->sigmaS = density * sigma_s;
        mRec->sigmaN = majorant - density * (sigma_a + sigma_s);
        mRec->tr = evalTransmittance(ray.origin, mRec->scatterPoint);
        // calculate pdf, i.e., sum of $1\n * \sigma_t^i e^{-\sigma_t^i * t}$.
        for (int i = 0; i < nSpectrumSamples; i++) {
            mRec->pdf[i] = majorant[i] * fm::exp(-majorant[i] * distance);
        }
        return true;
    } else {
        // sampled a point on object boundary (surface).
        mRec->marchLength = its.t;
        mRec->scatterPoint = ray.at(its.t);

        Point3d index = worldToIndex(ray.at(its.t));
        float density = sampleFromGrid(index, densityGrid);

        mRec->sigmaA = density * sigma_a;
        mRec->sigmaS = density * sigma_s;
        mRec->sigmaN = majorant - density * (sigma_a + sigma_s);
        mRec->tr = evalTransmittance(ray.origin, its.position);
        // calculate discrete probility (instead of continuous probability density), i.e., sum of $1\n * e^{-\sigma_t^i * t_max}$.
        for (int i = 0; i < nSpectrumSamples; i++) {
            mRec->pdf[i] = fm::exp(-majorant[i] * its.t);
        }
        return false;
    }
}

Spectrum NullScatteringMedium::evalTransmittance(Point3d from, Point3d dest) const {
    double dist = (dest - from).length();
    Spectrum tr;
    for (int i = 0; i < nSpectrumSamples; i++) {
        tr[i] = fm::exp(-majorant[i] * dist);
    }
    return tr;
}

MediumSampleEvent NullScatteringMedium::sampleEvent(Point3d pos,
                                                    Point2d sample) const {
    auto [x, y] = sample;
    int channelIndex = int(x * nSpectrumSamples);

    Point3d index = worldToIndex(pos);
    float density = sampleFromGrid(index, densityGrid);

    double a = density * sigma_a[channelIndex];
    double s = density * sigma_s[channelIndex];
    double t = majorant[channelIndex];

    if (y < a / t) {
        return MediumSampleEvent::a;
    } else if (y < (a + s) / t) {
        return MediumSampleEvent::s;
    } else {
        return MediumSampleEvent::n;
    }
}

Spectrum NullScatteringMedium::evalEmittance(Point3d pos) const {
    Point3d index = worldToIndex(pos);
    return sampleFromGrid(index, emissionGrid);
}

float NullScatteringMedium::sampleFromGrid(Point3d index,
                                           const nanovdb::FloatGrid *grid) const {
    if (grid == nullptr) {
        return .0f;
    }
    nanovdb::Vec3f indexf{float(index[0]), float(index[1]), float(index[2])};
    auto result = FloatGridSampler(grid->tree())(indexf);
    return result;
}

Spectrum NullScatteringMedium::sampleFromGrid(Point3d index,
                                              const nanovdb::Vec3fGrid *grid) const {
    if (grid == nullptr) {
        return Spectrum{.0};
    }
    nanovdb::Vec3f indexf{float(index[0]), float(index[1]), float(index[2])};
    auto result = Vec3fGridSampler(grid->tree())(indexf);
    return Spectrum(RGB3(result[0], result[1], result[2]));
}

Point3d NullScatteringMedium::worldToIndex(Point3d world) const {
    world = invTransformMatrix * world;
    auto index = densityGrid->worldToIndexF(nanovdb::Vec3f{(float)world[0], (float)world[1], (float)world[2]});
    return Point3d{index[0], index[1], index[2]};
}

Point3d NullScatteringMedium::indexToWorld(Point3d index) const {
    nanovdb::Vec3f idx{(float)index[0], (float)index[1], (float)index[2]};
    auto wrd = densityGrid->indexToWorldF(idx);
    Point3d world{wrd[0], wrd[1], wrd[2]};
    return transformMatrix * world;
}

Vec3d NullScatteringMedium::worldToIndexDir(Vec3d world) const {
    world = invTransformMatrix * world;
    auto index = densityGrid->worldToIndexDirF(nanovdb::Vec3f{(float)world[0], (float)world[1], (float)world[2]});
    return normalize(Vec3d{index[0], index[1], index[2]});
}

Vec3d NullScatteringMedium::indexToWorldDir(Vec3d index) const {
    nanovdb::Vec3f idx{(float)index[0], (float)index[1], (float)index[2]};
    auto wrd = densityGrid->indexToWorldDirF(idx);
    Vec3d world{wrd[0], wrd[1], wrd[2]};
    return normalize(transformMatrix * world);
}
