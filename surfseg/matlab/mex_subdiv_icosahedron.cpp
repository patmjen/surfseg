#include "mex.h"
#include "matrix.h"

#include "subdivided_icosahedron.h"
#include "matlab_util.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nsubdiv = nrhs > 0 ? getCastScalarChecked<int>(prhs[0], "nsubdiv") : 0;
	ensureOrError(nsubdiv >= 0, "nsubdiv must be nonnegative");

	SubdividedIcosahedron mesh(1.0f);
	mesh.subdivide(nsubdiv);
	size_t numVerts = mesh.vertices.size();
	size_t numFaces = mesh.faces.size();

	mxArray *vertices = mxCreateNumericMatrix(numVerts, 3, mxSINGLE_CLASS, mxREAL);
	mxArray *faces = mxCreateNumericMatrix(numFaces, 3, mxINT32_CLASS, mxREAL);
	float *vertData = static_cast<float *>(mxGetData(vertices));
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

	plhs[0] = faces;
	plhs[1] = vertices;

	// If requested, also compute and return vertex normals
	if (nlhs > 2) {
		mesh.computeVertexNormals();
		mxArray *normals = mxCreateNumericMatrix(numVerts, 3, mxSINGLE_CLASS, mxREAL);
		float *normData = static_cast<float *>(mxGetData(normals));
		for (const auto& v : mesh.vertices) {
			normData[v.self + 0 * numVerts] = v.normal[0];
			normData[v.self + 1 * numVerts] = v.normal[1];
			normData[v.self + 2 * numVerts] = v.normal[2];
		}
		plhs[2] = normals;
	}
}
