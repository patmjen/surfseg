#include <vector>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec3f.h>

#include "volume.h"
#include "surface_segment.h"
#include "matlab_util.h"

using namespace CGLA;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ensureOrError(nrhs == 10, "Must supply 11 inputs");
	Volume<float> cost = getVolumeChecked<float>(prhs[0], "Cost volume");
	Vec3f center = getVectorChecked<Vec3f>(prhs[1], "Center");
	float initR = getCastScalarChecked<float>(prhs[2], "initR");
	int initSubDiv = getCastScalarChecked<int>(prhs[3], "initSubDiv");
	int numSamples = getCastScalarChecked<int>(prhs[4], "numSamples");
	float sampleStep = getCastScalarChecked<float>(prhs[5], "sampleStep");
	int maxDiff = getCastScalarChecked<int>(prhs[6], "maxDiff");
	int minSep = getCastScalarChecked<int>(prhs[7], "minSep");
	int maxSep = getCastScalarChecked<int>(prhs[8], "maxSep");
	int k = getCastScalarChecked<int>(prhs[9], "k");
	CostType costType = static_cast<CostType>(getCastScalarChecked<int>(prhs[10], "costType"));

	// TODO: Extra input validation

	std::vector<Volume<float>> costs(k, cost);
	std::vector<ManifoldMesh> meshes = kSurfaceCutNoOverlap(costs, SubdividedIcosahedron(center, initR, initSubDiv),
		numSamples, sampleStep, maxDiff, minSep, maxSep, costType);

	mxArray *vertCell = mxCreateCellMatrix(1, k);
	mxArray *faceCell = mxCreateCellMatrix(1, k);

	for (int i = 0; i < k; ++i) {
		const ManifoldMesh& mesh = meshes[i];

		const size_t numVerts = mesh.vertices.size();
		const size_t numFaces = mesh.faces.size();
		mxArray *vertices = mxCreateNumericMatrix(numVerts, 3, mxSINGLE_CLASS, mxREAL);
		float *vertData = static_cast<float *>(mxGetData(vertices));
		mxArray *faces = mxCreateNumericMatrix(numFaces, 3, mxINT32_CLASS, mxREAL);
		int *faceData = static_cast<int *>(mxGetData(faces));

		for (const auto& v : mesh.vertices) {
			vertData[v.self + 0 * numVerts] = v.pos[0];
			vertData[v.self + 1 * numVerts] = v.pos[1];
			vertData[v.self + 2 * numVerts] = v.pos[2];
		}

		for (const auto& f : mesh.faces) {
			// Add one since MATLAB uses 1-indexing
			faceData[f.self + 0 * numFaces] = mesh.edges[f.edge].vert + 1;
			faceData[f.self + 1 * numFaces] = mesh.edges[mesh.next(f.edge)].vert + 1;
			faceData[f.self + 2 * numFaces] = mesh.edges[mesh.next(mesh.next(f.edge))].vert + 1;
		}

		mxSetCell(vertCell, i, vertices);
		mxSetCell(faceCell, i, faces);
	}

	plhs[0] = faceCell;
	plhs[1] = vertCell;
}
