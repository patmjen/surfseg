#include <vector>
#include <thread>
#include <memory>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec3f.h>

#include "volume.h"
#include "subdivided_icosahedron.h"
#include "surface_segment.h"
#include "matlab_util.h"

using namespace CGLA;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ensureOrError(9 <= nrhs && nrhs <= 11, "Must supply between 7 and 11 inputs");
	Volume<float> cost = getVolumeChecked<float>(prhs[0], "Cost volume");
	Volume<float> centers = getVolumeChecked<float>(prhs[1], "Centers");
	std::shared_ptr<float> initRs = getCastSharedPtrChecked<float>(prhs[2], "initR");
	int initSubDiv = getCastScalarChecked<int>(prhs[3], "initSubDiv");
	int numSamples = getCastScalarChecked<int>(prhs[4], "numSamples");
	float sampleStep = getCastScalarChecked<float>(prhs[5], "sampleStep");
	int maxDiff = getCastScalarChecked<int>(prhs[6], "maxDiff");
	const mxArray *connections = prhs[7];
	CostType costType = static_cast<CostType>(getCastScalarChecked<int>(prhs[8], "costType"));
	float smoothFactor = nrhs > 8 ? getCastScalarChecked<float>(prhs[9], "smoothFactor") : 0.0f;
	int smoothIter = nrhs > 9 ? getCastScalarChecked<int>(prhs[10], "smoothIter") : 0;

	// TODO: Extra input validation

	ensureOrError(centers.ny == 3 && isMatrix(prhs[1]), "Centers must be N x 3");
	ensureOrError(mxIsScalar(prhs[2]) || isVector(prhs[2]), "initR must be a scalar or vector");
	ensureOrError(mxIsSparse(connections), "Connections must be a sparse matrix");
	ensureOrError(!mxIsComplex(connections), "Connections must be real");

	int nr = mxGetNumberOfElements(prhs[2]);
	int nmesh = centers.nx;
	ensureOrError(nr == 1 || nr == nmesh, "Must provide a single initR or one for every center");
	int nthreads = std::min(getMaxCompThreads(), nmesh);
	int meshesPerThread = (nmesh + (nmesh % nthreads == 0 ? 0 : nthreads)) / nthreads;

	std::vector<SubdividedIcosahedron> meshes;
	std::vector<float *> vertData;
	std::vector<int *> faceData;

	// Make worker function to process subset
	auto workerFunc = [&](int begin, int end) {
		for (int i = begin; i < end; ++i) {
			// Extract needed planes
			size_t *ir = mxGetIr(connections);
			size_t *jc = mxGetJc(connections);
			size_t nconn = jc[i + 1] - jc[i]; // Number of nonzero entries in this column

			Vec3f center = meshes[i].center();
			float r = meshes[i].r();

			std::vector<Vec3f> planeNormals;
			std::vector<float> planeDists;
			planeNormals.reserve(nconn);
			planeDists.reserve(nconn);
			for (size_t ci = jc[i]; ci < jc[i+1]; ++ci) {
				size_t conn = ir[ci];
				Vec3f connCenter = meshes[conn].center();
				float connR = meshes[conn].r();

				Vec3f planeNorm = normalize(connCenter - center);
				float a = r / (r + meshes[conn].r());
				Vec3f planePos = (1.0f - a) * center + a * connCenter;
				float planeDist = dot(planeNorm, planePos);

				planeNormals.push_back(planeNorm);
				planeDists.push_back(planeDist);
			}

			ManifoldMesh mesh = surfaceCutPlaneSep(cost, std::move(meshes[i]), 
				numSamples, sampleStep, maxDiff, costType, planeNormals, planeDists);

			mesh.taubinSmooth(smoothFactor, 1.04f * smoothFactor, smoothIter);

			const size_t numVerts = mesh.vertices.size();
			const size_t numFaces = mesh.faces.size();

			for (const auto& v : mesh.vertices) {
				vertData[i][v.self + 0 * numVerts] = mesh.vpos(v)[0];
				vertData[i][v.self + 1 * numVerts] = mesh.vpos(v)[1];
				vertData[i][v.self + 2 * numVerts] = mesh.vpos(v)[2];
			}

			for (const auto& f : mesh.faces) {
				// Add one since MATLAB uses 1-indexing
				faceData[i][f.self + 0 * numFaces] = mesh.edges[f.edge].vert + 1;
				faceData[i][f.self + 1 * numFaces] = mesh.edges[mesh.next(f.edge)].vert + 1;
				faceData[i][f.self + 2 * numFaces] = mesh.edges[mesh.next(mesh.next(f.edge))].vert + 1;
			}
		}
	};

	// Allocate all needed memory
	mxArray *vertCell = mxCreateCellMatrix(1, nmesh);
	mxArray *faceCell = mxCreateCellMatrix(1, nmesh);
	meshes.reserve(nmesh);
	vertData.reserve(nmesh);
	faceData.reserve(nmesh);
	float *initR = initRs.get();
	for (int i = 0; i < nmesh; ++i) {
		Vec3f center(centers.at(i, 0, 0), centers.at(i, 1, 0), centers.at(i, 2, 0));
		meshes.push_back(SubdividedIcosahedron(center, *initR, initSubDiv));

		const size_t numVerts = meshes[i].vertices.size();
		const size_t numFaces = meshes[i].faces.size();
		mxArray *vertices = mxCreateNumericMatrix(numVerts, 3, mxSINGLE_CLASS, mxREAL);
		mxArray *faces = mxCreateNumericMatrix(numFaces, 3, mxINT32_CLASS, mxREAL);
		vertData.push_back(static_cast<float *>(mxGetData(vertices)));
		faceData.push_back(static_cast<int *>(mxGetData(faces)));

		mxSetCell(vertCell, i, vertices);
		mxSetCell(faceCell, i, faces);

		if (nr > 1) {
			initR++;
		}
	}

	// Compute all cuts
	if (nthreads == 1) {
		// Just run on this thread
		workerFunc(0, nmesh);
	} else {
		std::vector<std::thread> threads;
		// Split processing among threads
		for (int i = 0; i < nthreads; ++i) {
			int begin = i * meshesPerThread;
			int end = std::min((i + 1) * meshesPerThread, nmesh);
			threads.push_back(std::thread(workerFunc, begin, end));
		}
		// Wait for all threads to finish
		for (auto& th : threads) {
			if (th.joinable()) {
				th.join();
			}
		}
	}

	plhs[0] = faceCell;
	plhs[1] = vertCell;
}
