#include <vector>
#include <thread>
#include <memory>
#include <cmath>
#include <algorithm>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec4f.h>

#include "tet_mesh_4d.h"
#include "volume4d.h"
#include "surface_segment.h"
#include "matlab_util.h"

using namespace CGLA;

float simpleRegionCost(const Volume4d<float>& vol, Vec4f pos)
{
	float x = vol.interp(pos);
	return 10*x - (1 - x);
}

float simpleOnSurfaceCost(const Volume4d<float>& vol, Vec4f pos)
{
	// Use simple central difference for gradient
	Vec4f zero = Vec4f(0);

	Vec4f dxp = zero; dxp[0] = pos[0] <= vol.nx - 2 ? 1 : 0;
	Vec4f dxn = zero; dxp[0] = pos[0] >= 1 ? -1 : 0;

	Vec4f dyp = zero; dyp[1] = pos[1] <= vol.ny - 2 ? 1 : 0;
	Vec4f dyn = zero; dyp[1] = pos[1] >= 1 ? -1 : 0;

	Vec4f dzp = zero; dzp[2] = pos[2] <= vol.nz - 2 ? 1 : 0;
	Vec4f dzn = zero; dzp[2] = pos[2] >= 1 ? -1 : 0;

	Vec3f grad = 0.5 * Vec3f(
		vol.interp(pos + dxp) - vol.interp(pos + dxn),
		vol.interp(pos + dyp) - vol.interp(pos + dyn),
		vol.interp(pos + dzp) - vol.interp(pos + dzn)
	);
	return -length(grad);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ensureOrError(7 <= nrhs && nrhs <= 8, "Must supply between 7 and 8 inputs");
	ensureOrError(isSize(prhs[1], { -1, 4 }), "Tets must be an N x 4 array");
	ensureOrError(isSize(prhs[2], { -1, 4 }), "Vertices must be an M x 4 array");
	Volume4d<float> vol = getVolume4dChecked<float>(prhs[0], "Volume");
	Volume<int> tetVol = getCastVolumeChecked<int>(prhs[1], "Tets");
	Volume<float> vertVol = getCastVolumeChecked<float>(prhs[2], "Vertices");
	int numSamples = getCastScalarChecked<int>(prhs[3], "numSamples");
	float sampleStep = getCastScalarChecked<float>(prhs[4], "sampleStep");
	int maxDiff = getCastScalarChecked<int>(prhs[5], "maxDiff");
	CostType costType = static_cast<CostType>(getCastScalarChecked<int>(prhs[6], "costType"));
	std::vector<int> frozenVerts;
	if (nrhs > 7) {
		// TODO: This ends up doing unnecessary allocation, copying, and freeing. Check added overhead.
		Volume<int> frozenVertVol = getCastVolumeChecked<int>(prhs[7], "FrozenVertices");
		frozenVerts.assign(frozenVertVol.data.get(), frozenVertVol.data.get() + frozenVertVol.numElem());
		std::transform(frozenVerts.begin(), frozenVerts.end(), frozenVerts.begin(),
			[](auto x) { return x - 1; }); // Subtract 1 since MATLAB uses 1-indexing
	}

	size_t nverts = vertVol.nx;
	size_t ntets = tetVol.nx;

	std::vector<Vec4f> vertices;
	vertices.reserve(nverts);
	for (int i = 0; i < nverts; ++i) {
		vertices.push_back(Vec4f(vertVol.at(i, 0, 0), vertVol.at(i, 1, 0), vertVol.at(i, 2, 0), vertVol.at(i, 3, 0)));
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

	if (costType == CostType::ON_SURFACE) {
		mesh = surfaceCut4d(vol, mesh, numSamples, sampleStep, maxDiff, costType, simpleOnSurfaceCost,
			frozenVerts);
	} else {
		mesh = surfaceCut4d(vol, mesh, numSamples, sampleStep, maxDiff, costType, simpleRegionCost,
			frozenVerts);
	}

	mxArray *mxVerts = mxCreateNumericMatrix(nverts, 4, mxSINGLE_CLASS, mxREAL);
	mxArray *mxTets = mxCreateNumericMatrix(ntets, 4, mxINT32_CLASS, mxREAL);
	mxArray *mxNormals = mxCreateNumericMatrix(nverts, 4, mxSINGLE_CLASS, mxREAL);
	float *vertData = static_cast<float *>(mxGetData(mxVerts));
	int *tetData = static_cast<int *>(mxGetData(mxTets));
	float *normalData = static_cast<float *>(mxGetData(mxNormals));

	for (const auto& v : mesh.vertices) {
		vertData[v.self + 0 * nverts] = v.pos[0];
		vertData[v.self + 1 * nverts] = v.pos[1];
		vertData[v.self + 2 * nverts] = v.pos[2];
		vertData[v.self + 3 * nverts] = v.pos[3];

		normalData[v.self + 0 * nverts] = v.normal[0];
		normalData[v.self + 1 * nverts] = v.normal[1];
		normalData[v.self + 2 * nverts] = v.normal[2];
		normalData[v.self + 3 * nverts] = v.normal[3];
	}

	for (const auto& t : mesh.tets) {
		// Add one since MATLAB uses 1-indexing
		tetData[t.self + 0 * ntets] = t.verts[0] + 1;
		tetData[t.self + 1 * ntets] = t.verts[1] + 1;
		tetData[t.self + 2 * ntets] = t.verts[2] + 1;
		tetData[t.self + 3 * ntets] = t.verts[3] + 1;

		/*normalData[t.self + 0 * ntets] = t.normal[0];
		normalData[t.self + 1 * ntets] = t.normal[1];
		normalData[t.self + 2 * ntets] = t.normal[2];
		normalData[t.self + 3 * ntets] = t.normal[3];*/
	}

	plhs[0] = mxTets;
	plhs[1] = mxVerts;
	plhs[2] = mxNormals;
}