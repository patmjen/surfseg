#include <vector>
#include <unordered_map>
#include <array>
#include <utility>
#include <limits>
#include <cmath>

#include "mex.h"

#include <GEL/CGLA/Vec4i.h>
#include <GEL/CGLA/Vec4f.h>

#include "util.h"
#include "volume.h"
#include "tet_mesh_4d.h"
#include "matlab_util.h"

using namespace CGLA;

constexpr float TOL = 100 * std::numeric_limits<float>::epsilon();
constexpr double DBL_NAN = std::numeric_limits<double>::quiet_NaN();

inline int robustSign(float x, float tol)
{
    return -static_cast<int>(x < -tol) + static_cast<int>(x > tol);
}

inline Vec4f getVert(const Volume<float>& vertVol, int vertIdx)
{
    return Vec4f(
        vertVol.at(vertIdx, 0, 0),
        vertVol.at(vertIdx, 1, 0),
        vertVol.at(vertIdx, 2, 0),
        vertVol.at(vertIdx, 3, 0)
    );
}

std::array<int, 4> sortIdx4(int v1, int v2, int v3, int v4)
{
    std::array<int, 4> vals = { v1, v2, v3, v4 };
    std::array<int, 4> idxs = { 0, 1, 2, 3 };
    // Just use bubble sort, since we need so few iterations
    for (int si = 0; si < 3; ++si) {
        for (int i = 0; i < 4 - 1 - si; ++i) {
            if (vals[i] > vals[i + 1]) {
                std::swap(vals[i], vals[i + 1]);
                std::swap(idxs[i], idxs[i + 1]);
            }
        }
    }
    return idxs;
}


std::pair<int, int> makeKey(int vk1, int vk2)
{
    return vk1 <= vk2 ? std::make_pair(vk1, vk2) : std::make_pair(vk2, vk1);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ensureOrError(3 <= nrhs && nrhs <= 5, "Must supply between 3 and 5 inputs");
	ensureOrError(isSize(prhs[0], { -1, 4 }), "Tets must be an N x 4 array");
	ensureOrError(isSize(prhs[1], { -1, 4 }), "Vertices must be an M x 4 array");
    const Volume<int> tetVol = getCastVolumeChecked<int>(prhs[0], "Tets");
	const Volume<float> vertVol = getCastVolumeChecked<float>(prhs[1], "Vertices");
    const float time = getCastScalarChecked<float>(prhs[2], "time");
    const Vec4f normal = nrhs < 4 ? Vec4f(0, 0, 0, 1) : getCastVectorChecked<Vec4f>(prhs[3], "Normal");
    const bool triMesh = nrhs < 5 ? false : getCastScalarChecked<bool>(prhs[4], "triMesh");

    const bool saveIntTets = nlhs > 2;

    std::unordered_map<std::pair<int, int>, int> edgeToPoint;
    std::vector<Vec4f> vertices;
    std::vector<Vec4i> faces; // For triangle meshes, we just ignore the last entry
    std::vector<int> intTets;
    std::vector<int> faceTetIdxs; // Index of tet. each face "came from"

    auto addVertex = [&](int vk, Vec4f v)
    {
        auto key = std::make_pair(vk, vk); // No need to check order, so just make the pair
        if (edgeToPoint.find(key) == edgeToPoint.end()) {
            // Key not seen before, so add the new vertex
            vertices.push_back(v);
            edgeToPoint[key] = vertices.size() - 1;
        }
        return edgeToPoint[key];
    };

    auto addIntersection = [&](int vk1, int vk2, Vec4f v1, Vec4f v2, int s1, int s2, float d1, float d2)
    {
        assert(s1 != s2); // Both vertices may **not** be on the same side

        if (s1 == 0) {
            // First vertex is in hyperplane
            return addVertex(vk1, v1);
        } else if (s2 == 0) {
            // Second vertex is in hyperplane
            return addVertex(vk2, v2);
        } else {
            // No vertex is in the hyperplane, so compute intersection
            auto key = makeKey(vk1, vk2);
            if (edgeToPoint.find(key) == edgeToPoint.end()) {
                // Pair not seen before, so add the new vertex
                float a = d1 / (d1 - d2);
                Vec4f p = (1 - a) * v1 + a * v2;
                vertices.push_back(p);
                edgeToPoint[key] = vertices.size() - 1;
            }
            return edgeToPoint[key];
        }
    };

    // Find all tets which intersects the hyperplane, and store their intersection
    for (int ti = 0; ti < tetVol.nx; ++ti) {
        // Get tet vertices (subtract 1 since MATLAB uses 1-indexing)
        int vk0 = tetVol.at(ti, 0, 0) - 1;
        int vk1 = tetVol.at(ti, 1, 0) - 1;
        int vk2 = tetVol.at(ti, 2, 0) - 1;
        int vk3 = tetVol.at(ti, 3, 0) - 1;
        Vec4f v0 = getVert(vertVol, vk0);
        Vec4f v1 = getVert(vertVol, vk1);
        Vec4f v2 = getVert(vertVol, vk2);
        Vec4f v3 = getVert(vertVol, vk3);

        // Check which vertices lie above and below the hyperplane
        float d0 = dot(v0, normal) - time;
        float d1 = dot(v1, normal) - time;
        float d2 = dot(v2, normal) - time;
        float d3 = dot(v3, normal) - time;
        int s0 = robustSign(d0, TOL);
        int s1 = robustSign(d1, TOL);
        int s2 = robustSign(d2, TOL);
        int s3 = robustSign(d3, TOL);
        int npos = (s0 > 0) + (s1 > 0) + (s2 > 0) + (s3 > 0);
        int nneg = (s0 < 0) + (s1 < 0) + (s2 < 0) + (s3 < 0);
        int nzer = 4 - npos - nneg;

        if (npos == 4 || nneg == 4) {
            // No intersection, move on to next
            continue;
        }
        intTets.push_back(ti); // Mark tet as intersecting
        if (nzer == 4) {
            // Intersection is whole tet
            int nvk0 = addVertex(vk0, v0);
            int nvk1 = addVertex(vk1, v1);
            int nvk2 = addVertex(vk2, v2);
            int nvk3 = addVertex(vk3, v3);
            faces.push_back(Vec4i(nvk0, nvk1, nvk2, -1));
            faces.push_back(Vec4i(nvk0, nvk1, nvk3, -1));
            faces.push_back(Vec4i(nvk0, nvk2, nvk3, -1));
            faces.push_back(Vec4i(nvk1, nvk2, nvk3, -1));

            faceTetIdxs.push_back(ti);
            faceTetIdxs.push_back(ti);
            faceTetIdxs.push_back(ti);
            faceTetIdxs.push_back(ti);
        } else if (npos == 2 && nneg == 2) {
            // Intersection is a quad
            const auto idxs = sortIdx4(s0, s1, s2, s3);

            const std::array<int, 4> tetKeys = { vk0, vk1, vk2, vk3 };
            const std::array<Vec4f, 4> tetVerts = { v0, v1, v2, v3 };
            const std::array<int, 4> tetSigns = { s0, s1, s2, s3 };
            const std::array<float, 4> tetDists = { d0, d1, d2, d3 };
            auto addTetIntersection = [&](int i1, int i2)
            {
                return addIntersection(
                    tetKeys[i1], tetKeys[i2],
                    tetVerts[i1], tetVerts[i2],
                    tetSigns[i1], tetSigns[i2],
                    tetDists[i1], tetDists[i2]
                );
            };

            int nvk0 = addTetIntersection(idxs[0], idxs[2]);
            int nvk1 = addTetIntersection(idxs[0], idxs[3]);
            int nvk2 = addTetIntersection(idxs[1], idxs[2]);
            int nvk3 = addTetIntersection(idxs[1], idxs[3]);
            Vec4f nv0 = vertices[nvk0];
            Vec4f nv1 = vertices[nvk1];
            Vec4f nv2 = vertices[nvk2];
            Vec4f nv3 = vertices[nvk3];
            Vec4f center = 0.25 * (nv0 + nv1 + nv2 + nv3);

            // For quads, we need to make sure the vertices are in clockwise order for rendering
            Vec4f nrm = cross4(nv0 - center, nv1 - center, normal);
            if (dot(nrm, cross4(nv2 - center, nv3 - center, normal)) < 0) {
                // Last two vertices not in clockwise order, so swap them
                std::swap(nvk2, nvk3);
            }

            if (triMesh) {
                faces.push_back(Vec4i(nvk0, nvk1, nvk2, -1));
                faces.push_back(Vec4i(nvk0, nvk2, nvk3, -1));

                faceTetIdxs.push_back(ti);
                faceTetIdxs.push_back(ti);
            } else {
                faces.push_back(Vec4i(nvk0, nvk1, nvk2, nvk3));

                faceTetIdxs.push_back(ti);
            }
        } else if (nzer == 1) {
            // Intersection is a triangular tet face
            const auto idxs = sortIdx4(abs(s0), abs(s1), abs(s2), abs(s3)); // We don't care about above or below

            const std::array<int, 4> tetKeys = { vk0, vk1, vk2, vk3 };
            const std::array<Vec4f, 4> tetVerts = { v0, v1, v2, v3 };
            std::array<int, 3> newKeys;
            for (int i = 0; i < 3; ++i) {
                newKeys[i] = addVertex(tetKeys[idxs[i]], tetVerts[idxs[i]]);
            }
            faces.push_back(Vec4i(newKeys[0], newKeys[1], newKeys[2], -1));

            faceTetIdxs.push_back(ti);
        } else if (npos == 1 || nneg == 1) {
            // Intersection is a triangle
            const auto idxs = sortIdx4(s0, s1, s2, s3);
            int idxA, idxB1, idxB2, idxB3;
            if (npos == 1) {
                idxA = idxs[3];
                idxB1 = idxs[0];
                idxB2 = idxs[1];
                idxB3 = idxs[2];
            } else {
                idxA = idxs[0];
                idxB1 = idxs[1];
                idxB2 = idxs[2];
                idxB3 = idxs[3];
            }

            const std::array<int, 4> tetKeys = { vk0, vk1, vk2, vk3 };
            const std::array<Vec4f, 4> tetVerts = { v0, v1, v2, v3 };
            const std::array<int, 4> tetSigns = { s0, s1, s2, s3 };
            const std::array<float, 4> tetDists = { d0, d1, d2, d3 };
            auto addTetIntersection = [&](int i1, int i2)
            {
                return addIntersection(
                    tetKeys[i1], tetKeys[i2],
                    tetVerts[i1], tetVerts[i2],
                    tetSigns[i1], tetSigns[i2],
                    tetDists[i1], tetDists[i2]
                );
            };

            int nvk0 = addTetIntersection(idxA, idxB1);
            int nvk1 = addTetIntersection(idxA, idxB2);
            int nvk2 = addTetIntersection(idxA, idxB3);

            faces.push_back(Vec4i(nvk0, nvk1, nvk2, -1));

            faceTetIdxs.push_back(ti);
        } else if (nzer == 2 && (npos == 2 || nneg == 2)) {
            // Intersection is line segment
            const auto idxs = sortIdx4(abs(s0), abs(s1), abs(s2), abs(s3)); // We don't care about above or below

            const std::array<int, 4> tetKeys = { vk0, vk1, vk2, vk3 };
            const std::array<Vec4f, 4> tetVerts = { v0, v1, v2, v3 };

            int nvk0 = addVertex(tetKeys[idxs[0]], tetVerts[idxs[0]]);
            int nvk1 = addVertex(tetKeys[idxs[1]], tetVerts[idxs[1]]);

            faces.push_back(Vec4i(nvk0, nvk1, -1, -1));

            faceTetIdxs.push_back(ti);
        } else if (nzer == 1 && (npos == 3 || nneg == 3)) {
            // Intersection is a point
            const auto idxs = sortIdx4(abs(s0), abs(s1), abs(s2), abs(s3)); // We don't care about above or below

            const std::array<int, 4> tetKeys = { vk0, vk1, vk2, vk3 };
            const std::array<Vec4f, 4> tetVerts = { v0, v1, v2, v3 };

            addVertex(tetKeys[idxs[0]], tetVerts[idxs[0]]);
        } else {
            ensureOrError(false, "Invalid tet. configuration"); // We should never hit this!
        }
    }

    // Put results into MATLAB arrays
    size_t numFaces = faces.size();
	mxArray *mxFaces = mxCreateNumericMatrix(numFaces, triMesh ? 3 : 4, mxDOUBLE_CLASS, mxREAL);
    double *faceData = static_cast<double *>(mxGetData(mxFaces));
    for (int i = 0; i < numFaces; ++i) {
        Vec4i f = faces[i];
        // Add 1 since MATLAB uses 1-indexing
        faceData[i + 0 * numFaces] = f[0] == -1 ? DBL_NAN : f[0] + 1;
        faceData[i + 1 * numFaces] = f[1] == -1 ? DBL_NAN : f[1] + 1;
        faceData[i + 2 * numFaces] = f[2] == -1 ? DBL_NAN : f[2] + 1;
        if (!triMesh) {
            faceData[i + 3 * numFaces] = f[3] == -1 ? DBL_NAN : f[3] + 1;
        }
    }
	plhs[0] = mxFaces;

    size_t numVerts = vertices.size();
    mxArray *mxVerts = mxCreateNumericMatrix(numVerts, 4, mxDOUBLE_CLASS, mxREAL);
    double *vertData = static_cast<double *>(mxGetData(mxVerts));
    for (int i = 0; i < numVerts; ++i) {
        Vec4f v = vertices[i];
        vertData[i + 0 * numVerts] = v[0];
        vertData[i + 1 * numVerts] = v[1];
        vertData[i + 2 * numVerts] = v[2];
        vertData[i + 3 * numVerts] = v[3];
    }
	plhs[1] = mxVerts;

    if (saveIntTets) {
        size_t numIntTets = intTets.size();
        mxArray *mxIntTets = mxCreateNumericMatrix(numIntTets, 1, mxDOUBLE_CLASS, mxREAL);
        double *intTetData = static_cast<double *>(mxGetData(mxIntTets));
        for (int i = 0; i < numIntTets; ++i) {
            // Add 1 since MATLAB uses 1-indexing
            intTetData[i] = intTets[i] + 1;
        }
        plhs[2] = mxIntTets;

        mxArray *mxFaceTetIdxs = mxCreateNumericMatrix(numFaces, 1, mxDOUBLE_CLASS, mxREAL);
        double *faceTetIdxsData = static_cast<double *>(mxGetData(mxFaceTetIdxs));
        for (int i = 0; i < numFaces; ++i) {
            // Add 1 since MATLAB uses 1-indexing
            faceTetIdxsData[i] = faceTetIdxs[i] + 1;
        }
        plhs[3] = mxFaceTetIdxs;
    }
}