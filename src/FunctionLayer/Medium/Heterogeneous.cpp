#include "Heterogeneous.h"
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>
#include <FunctionLayer/Sampler/Independent.h>
#include <ResourceLayer/File/FileUtils.h>

class RegularTracker {
public:
    RegularTracker() = delete;

    RegularTracker(const Vec3i &min,
                   const Vec3i &max,
                   Point3d origin,
                   Vec3d direction,
                   float _tmax,
                   float _voxelSize,
                   float *t_world = nullptr) {
        voxelSize = _voxelSize * 0.03;
        tmin = .0f;
        tmax = _tmax / voxelSize;

        {
            // Check if the ray overlap the voxelGrid
            float minBound = -FLT_MAX, maxBound = FLT_MAX;
            for (int axis = 0; axis < 3; ++axis) {
                float t1 = (min[axis] - origin[axis]) / direction[axis],
                      t2 = (max[axis] - origin[axis]) / direction[axis];
                if (t1 > t2) std::swap(t1, t2);
                minBound = std::max(minBound, t1);
                maxBound = std::min(maxBound, t2);

                if ((minBound > maxBound) || maxBound < 0) terminate = true;
            }

            if (!terminate) {
                // Set the ray origin at the bound of the grid
                float t = minBound > 0 ? minBound : 0;
                origin += direction * t;
                tmax -= t;
                if (t_world) *t_world = t * voxelSize;
            }
        }

        for (int axis = 0; axis < 3; ++axis) {

            voxel[axis] = clamp((int)std::floor(origin[axis]), min[axis], max[axis]);
            deltaT[axis] = 1.f / std::abs(direction[axis]);

            if (direction[axis] == -.0f) direction[axis] = .0f;
            if (direction[axis] >= 0) {
                nextCrossingT[axis] = (voxel[axis] + 1.f - origin[axis]) / direction[axis];
                step[axis] = 1;
                voxelLimit[axis] = max[axis] + 1;
            } else {
                nextCrossingT[axis] = (voxel[axis] - origin[axis]) / direction[axis];
                step[axis] = -1;
                voxelLimit[axis] = min[axis] - 1;
            }
        }
    }

    bool track(Vec3i &index, float *dt) {
        if (terminate) return false;

        int stepAxis = -1;

        if (nextCrossingT[0] < nextCrossingT[1] && nextCrossingT[0] < nextCrossingT[2])
            stepAxis = 0;
        else if (nextCrossingT[1] < nextCrossingT[2])
            stepAxis = 1;
        else
            stepAxis = 2;

        if (nextCrossingT[stepAxis] > tmax) {
            *dt = (tmax - tmin) * voxelSize;
            tmin = tmax;
            terminate = true;
            //     std::cout << "True";
        } else {
            *dt = (nextCrossingT[stepAxis] - tmin) * voxelSize;
            tmin = nextCrossingT[stepAxis];
        }

        index[0] = voxel[0];
        index[1] = voxel[1];
        index[2] = voxel[2];

        voxel[stepAxis] += step[stepAxis];

        if (voxel[stepAxis] == voxelLimit[stepAxis])
            terminate = true;

        nextCrossingT[stepAxis] += deltaT[stepAxis];
        return true;
    }

public:
    bool terminate = false;
    float voxelSize;
    float tmin, tmax;
    float nextCrossingT[3], deltaT[3];
    int step[3], voxelLimit[3], voxel[3];
};

using FloatGridSampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 0, false>;
using Vec3fGridSampler = nanovdb::SampleFromVoxels<nanovdb::Vec3fGrid::TreeType, 0, false>;

Point3d HeterogeneousMedium::worldToIndex(Point3d world) const {
    world = invTransformMatrix * world;
    auto index = densityGrid->worldToIndexF(nanovdb::Vec3f{(float)world[0], (float)world[1], (float)world[2]});
    return Point3d{index[0], index[1], index[2]};
}

Point3d HeterogeneousMedium::indexToWorld(Point3d index) const {
    nanovdb::Vec3f idx{(float)index[0], (float)index[1], (float)index[2]};
    auto wrd = densityGrid->indexToWorldF(idx);
    Point3d world{wrd[0], wrd[1], wrd[2]};
    return transformMatrix * world;
}

Vec3d HeterogeneousMedium::worldToIndexDir(Vec3d world) const {
    world = invTransformMatrix * world;
    auto index = densityGrid->worldToIndexDirF(nanovdb::Vec3f{(float)world[0], (float)world[1], (float)world[2]});
    return normalize(Vec3d{index[0], index[1], index[2]});
}

Vec3d HeterogeneousMedium::indexToWorldDir(Vec3d index) const {
    nanovdb::Vec3f idx{(float)index[0], (float)index[1], (float)index[2]};
    auto wrd = densityGrid->indexToWorldDirF(idx);
    Vec3d world{wrd[0], wrd[1], wrd[2]};
    return normalize(transformMatrix * world);
}

HeterogeneousMedium::HeterogeneousMedium(std::string _gridFilePath,
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
}

bool HeterogeneousMedium::sampleDistance(MediumSampleRecord *mRec, const Ray &ray, const Intersection &its, Point2d sample) const {

    auto [x, y] = sample;

    // * randomly pick a channel/frequency and sample a distance.
    int channelIndex = int(x * nSpectrumSamples);

    Point3d origin = worldToIndex(ray.origin);
    Vec3d direction = worldToIndexDir(ray.direction);

    Vec3i index;
    float thick = -fm::log(1 - y);
    float dt;
    float sum = 0.f;
    float t_world = 0.f;

    RegularTracker rt(minIndex, maxIndex, origin, direction, its.t, voxelSize, &t_world);

    while (rt.track(index, &dt)) {
        Point3d voxel_loc(index[0] + .5f, index[1] + .5f, index[2] + .5f);
        float density = sampleFromGrid(voxel_loc, densityGrid) * (sigma_a + sigma_s)[channelIndex];
        float delta = density * dt;

        if (sum + delta >= thick) {
            // Sample the collision point

            dt = (thick - sum) / density;
            t_world += dt;
            mRec->scatterPoint = ray.at(t_world);

            mRec->marchLength = t_world;
            mRec->sigmaA = sigma_a * density;
            mRec->sigmaS = sigma_s * density;
            mRec->tr = Spectrum(fm::exp(-thick));
            mRec->pdf = mRec->tr * (sigma_a + sigma_s) * density;// TODO

            return true;
        }
        sum += delta;
        t_world += dt;
    }

    mRec->marchLength = its.t;
    mRec->tr = Spectrum(fm::exp(-sum));
    mRec->pdf = mRec->tr;// TODO
    return false;
}

Spectrum HeterogeneousMedium::evalTransmittance(Point3d from, Point3d dest) const {
    ////Point3d origin = from;
    ////Vec3d direction = normalize(dest - from);
    ////
    ////auto o_grid = densityGrid->worldToIndexF(nanovdb::Vec3f(origin[0], origin[1], origin[2])),
    ////     d_grid = densityGrid->worldToIndexDirF(nanovdb::Vec3f(direction[0], direction[1], direction[2]));
    ////
    ////origin = Point3d{o_grid[0], o_grid[1], o_grid[2]};
    ////direction = normalize(Vec3d(d_grid[0], d_grid[1], d_grid[2]));

    Vec3d direction = normalize(dest - from);
    double t = (dest - from).length();
    Point3d origin = worldToIndex(from);

    RegularTracker rt(minIndex, maxIndex, origin, direction, t, voxelSize);

    Vec3i index;
    float dt;
    float thick = .0f;
    while (rt.track(index, &dt)) {
        Point3d voxel_loc(index[0] + .5, index[1] + .5, index[2] + .5);
        float density = sampleFromGrid(voxel_loc, densityGrid);
        thick += density * dt;
    }
    return Spectrum(RGB3(fm::exp(-thick * (sigma_a + sigma_s)[0]),
                         fm::exp(-thick * (sigma_a + sigma_s)[1]),
                         fm::exp(-thick * (sigma_a + sigma_s)[2])));
}

float HeterogeneousMedium::sampleFromGrid(Point3d index,
                                          const nanovdb::FloatGrid *grid) const {
    if (grid == nullptr) {
        return .0f;
    }
    nanovdb::Vec3f indexf{float(index[0]), float(index[1]), float(index[2])};
    auto result = FloatGridSampler(grid->tree())(indexf);
    return result;
}

Spectrum HeterogeneousMedium::sampleFromGrid(Point3d index,
                                             const nanovdb::Vec3fGrid *grid) const {
    if (grid == nullptr) {
        return Spectrum{.0};
    }
    nanovdb::Vec3f indexf{float(index[0]), float(index[1]), float(index[2])};
    auto result = Vec3fGridSampler(grid->tree())(indexf);
    return Spectrum(RGB3(result[0], result[1], result[2]));
}

MediumSampleEvent HeterogeneousMedium::sampleEvent(Point3d pos,
                                                   Point2d sample) const {
    auto [x, y] = sample;
    int channelIndex = int(x * nSpectrumSamples);

    Point3d index = worldToIndex(pos);
    float density = sampleFromGrid(index, densityGrid);

    float a = density * sigma_a[channelIndex];
    float s = density * sigma_s[channelIndex];

    if (y < a / (a + s)) {
        return MediumSampleEvent::a;
    } else {
        return MediumSampleEvent::s;
    }
}

Spectrum HeterogeneousMedium::evalEmittance(Point3d pos) const {
    Point3d index = worldToIndex(pos);
    return sampleFromGrid(index, emissionGrid);
}
