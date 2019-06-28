#ifndef DSC_H__
#define DSC_H__

#include <vector>

#include <GEL/CGLA/Vec3f.h>

#include "simplex_mesh.h"

using namespace CGLA;

/** Deformable Simplicial Complex */
class DSC : public SimplexMesh {
public:
    explicit DSC() = default;

    float tetQuality(TetKey tk) const;
    std::vector<float> allTetQuality() const;

    void smartLaplacianSmooth(bool interfaceFlagsSet = true);
    void smartLaplacianSmooth(std::vector<float>& tetQualities,
        bool interfaceFlagsSet = true);
};

#endif // DSC_H__
