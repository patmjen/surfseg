#include <vector>
#include <unordered_set>
#include <utility>
#include <array>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Mat3x3d.h>

#include "tet_mesh_4d.h"
#include "util.h"
#include "matlab_util.h"

using namespace CGLA;

using VertKey = TetMesh4d::VertKey;
using Tetrahedron = TetMesh4d::Tetrahedron;

std::pair<std::array<VertKey, 3>, int> tetTriVertOrder(const Tetrahedron & tet, 
    const std::array<VertKey, 3>& triVerts)
{
    int oi = 0;
    int removedVert = TetMesh4d::INVALID_VERTEX;
    int vi = 0;
    std::array<VertKey, 3> out;
    for (auto vk : tet.verts) {
        int i;
        for (i = 0; i < 3; ++i) {
            if (vk == triVerts[i]) {
                out[oi] = vk;
                ++oi;
                break;
            }
        }
        if (i == 3) {
            removedVert = vi;
        }
        ++vi;
    }
    return std::make_pair(out, removedVert);
}

inline std::pair<int, int> makePairKey(int x, int y)
{
    return x < y ? std::make_pair(x, y) : std::make_pair(y, x);
}

inline bool isEven(int x)
{
    return x % 2 == 0;
}

// Compute (-1)^p
inline double negOnePow(int p) {
    return isEven(p) ? 1 : -1;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    ensureArgCount(nrhs, 2);
	Volume<int> tetVol = getCastVolumeChecked<int>(prhs[0], "Tets");
	Volume<float> vertVol = getCastVolumeChecked<float>(prhs[1], "Vertices");

    size_t nverts = vertVol.nx;
    size_t ntets = tetVol.nx;

    std::vector<Vec4f> vertices;
    vertices.reserve(nverts);
    for (int i = 0; i < nverts; ++i) {
        vertices.push_back(
            Vec4f(vertVol.at(i, 0, 0), vertVol.at(i, 1, 0), vertVol.at(i, 2, 0), vertVol.at(i, 3, 0))
        );
    }
    std::vector<std::array<int, 4>> tets;
    tets.reserve(ntets);
    for (int i = 0; i < ntets; ++i) {
        // Subtract 1 since MATLAB uses 1-indexing
        tets.push_back({
            tetVol.at(i, 0, 0) - 1,
            tetVol.at(i, 1, 0) - 1,
            tetVol.at(i, 2, 0) - 1,
            tetVol.at(i, 3, 0) - 1
        });
    }
    TetMesh4d mesh(vertices, tets);

    std::unordered_set<std::pair<int, int>> inconsistentTets;
    for (const auto& tet : mesh.tets) {
        for (auto tk : tet.tris) {
            const auto& tri = mesh.triangles[tk];
            const auto adjTetKey = tri.t0 == tet.self ? tri.t1 : tri.t0;
            if (adjTetKey == TetMesh4d::INVALID_TETRAHEDRON) {
                // No adjacent tet. for this triangle, so just move on
                continue;
            }
            const auto& adjTet = mesh.tets[adjTetKey];
            auto triVerts = mesh.getTriVerts(tk);
            std::array<VertKey, 3> tetOrder, adjOrder; 
            int tetRemoved, adjRemoved;
            std::tie(tetOrder, tetRemoved) = tetTriVertOrder(tet, triVerts);
            std::tie(adjOrder, adjRemoved) = tetTriVertOrder(adjTet, triVerts);

            Mat3x3d M(0);
            for (int i = 0; i < 3; ++i) {
                int j;
                for (j = 0; tetOrder[i] != adjOrder[j]; ++j) {}
                M[i][j] = 1;
            }
            if (determinant(M)*negOnePow(tetRemoved)*negOnePow(adjRemoved) != -1) {
                if (nlhs > 1) {
                    inconsistentTets.insert(makePairKey(tet.self + 1, adjTet.self + 1));
                }
            }
        }
    }

    plhs[0] = mxCreateDoubleScalar(static_cast<double>(inconsistentTets.size()));
    if (nlhs > 1) {
        mxArray *mxInConTets = mxCreateNumericMatrix(inconsistentTets.size(), 2, mxINT32_CLASS, mxREAL);
        int *ictData = static_cast<int *>(mxGetData(mxInConTets));
        int i = 0;
        for (const auto& ict : inconsistentTets) {
            std::tie(ictData[i], ictData[i + inconsistentTets.size()]) = ict;
            ++i;
        }
        plhs[1] = mxInConTets;
    }
}