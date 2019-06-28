#include <vector>
#include <cstring>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec3f.h>

#include "volume.h"
#include "manifold_mesh.h"
#include "matlab_util.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ensureOrError(nrhs == 2, "Must supply 2 inputs");
    ensureOrError(isSize(prhs[0], { -1, 3 }), "Faces must be an N x 3 array");
    ensureOrError(isSize(prhs[1], { -1, 3 }), "Vertices must be an M x 3 array");
    Volume<int> faceVol = getVolumeChecked<int>(prhs[0], "Faces");
    Volume<float> vertVol = getVolumeChecked<float>(prhs[1], "Vertices");

    size_t nverts = vertVol.nx;
    size_t nfaces = faceVol.nx;

    std::vector<Vec3f> vertices;
    vertices.reserve(nverts);
    for (int i = 0; i < nverts; ++i) {
        vertices.push_back(Vec3f(vertVol.at(i, 0, 0), vertVol.at(i, 1, 0), vertVol.at(i, 2, 0)));
    }
    std::vector<std::vector<int>> faces;
    faces.reserve(nfaces);
    for (int i = 0; i < nfaces; ++i) {
        faces.push_back({ faceVol.at(i, 0, 0), faceVol.at(i, 1, 0), faceVol.at(i, 2, 0) });
    }

    ManifoldMesh mesh;
    mesh.build(vertices, faces);

    mxArray *mxGaussCurvs = mxCreateNumericMatrix(1, vertVol.nx, mxSINGLE_CLASS, mxREAL);
    std::vector<float> gaussCurvs = mesh.gaussCurvatures();
    std::memcpy(mxGetData(mxGaussCurvs), gaussCurvs.data(), nverts * sizeof(float));

    plhs[0] = mxGaussCurvs;
}