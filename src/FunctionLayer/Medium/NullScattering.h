#pragma once

#include <nanovdb/util/GridHandle.h>
#include "Medium.h"
#include <CoreLayer/Geometry/Matrix.h>
using BufferT = nanovdb::HostBuffer;

class NullScatteringMedium : public Medium {
public:
    NullScatteringMedium(std::string gridFilePath,
                         std::shared_ptr<PhaseFunction> phase,
                         TransformMatrix3D _transform,
                         Spectrum _sigma_a,
                         Spectrum _sigma_s);

    virtual bool sampleDistance(MediumSampleRecord *mRec,
                                const Ray &ray,
                                const Intersection &its,
                                Point2d sample) const override;

    virtual Spectrum evalTransmittance(Point3d from,
                                       Point3d dest) const override;

    virtual MediumSampleEvent sampleEvent(Point3d pos,
                                          Point2d sample) const override;

    virtual Spectrum evalEmittance(Point3d pos) const override;

protected:
    float sampleFromGrid(Point3d index,
                         const nanovdb::FloatGrid *grid) const;

    Spectrum sampleFromGrid(Point3d index,
                            const nanovdb::Vec3fGrid *grid) const;

    Point3d worldToIndex(Point3d world) const;

    Point3d indexToWorld(Point3d index) const;

    Vec3d worldToIndexDir(Vec3d world) const;

    Vec3d indexToWorldDir(Vec3d index) const;

private:
    mutable TransformMatrix3D transformMatrix, invTransformMatrix;
    Vec3i minIndex, maxIndex;
    float voxelSize;

    nanovdb::GridHandle<BufferT> densityGridBuffer, emissionGridBuffer;
    const nanovdb::FloatGrid *densityGrid = nullptr;
    const nanovdb::Vec3fGrid *emissionGrid = nullptr;

    Spectrum sigma_a;
    Spectrum sigma_s;
    Spectrum majorant;// 由于全局使用同一个majorant值，因此在介质密度极不均衡时会产生较大的方差。
};
