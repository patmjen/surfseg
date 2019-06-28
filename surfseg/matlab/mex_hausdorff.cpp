#include <cmath>

#include "mex.h"
#include "matrix.h"

#include "compute_error.h"
#include "geomutils.h"
#undef max
#undef min

#include "matlab_util.h"
#include "volume.h"
#include "util.h"

inline vertex_t vmax(const vertex_t& v1, const vertex_t& v2)
{
	vertex_t out;
	out.x = std::fmax(v1.x, v2.x);
	out.y = std::fmax(v1.y, v2.y);
	out.z = std::fmax(v1.z, v2.z);
	return out;
}

inline vertex_t vmin(const vertex_t& v1, const vertex_t& v2)
{
	vertex_t out;
	out.x = std::fmin(v1.x, v2.x);
	out.y = std::fmin(v1.y, v2.y);
	out.z = std::fmin(v1.z, v2.z);
	return out;
}

model_error initModelError()
{
	model_error me;
	me.min_error = 0;
	me.max_error = 0;
	me.mean_error = 0;
	me.mesh = nullptr;
	me.n_samples = 0;
	me.fe = nullptr;
	me.verror = nullptr;
	me.info = nullptr;

	return me;
}

dist_surf_surf_stats initDistSurfSurfStats()
{
	size3d grid_size;
	grid_size.x = 0;
	grid_size.y = 0;
	grid_size.z = 0;

	dist_surf_surf_stats s;
	s.st_m1_area = 0;
	s.m1_area = 0;
	s.m2_area = 0;
	s.min_dist = 0;
	s.max_dist = 0;
	s.mean_dist = 0;
	s.rms_dist = 0;
	s.cell_sz = 0;
	s.n_t_p_nec = 0;
	s.m1_samples = 0;
	s.grid_sz = grid_size;
	s.n_ne_cells = 0;

	return s;
}

model buildModel(const Volume<int> faces, const Volume<float> verts)
{
	// Init fields and alloc memory
	int nfaces = faces.nx;
	int nverts = verts.nx;
	model out;
	out.num_faces = nfaces;
	out.num_vert = nverts;
	out.faces = static_cast<face_t *>(mxMalloc(nfaces * sizeof(face_t)));
	out.vertices = static_cast<vertex_t *>(mxMalloc(nverts * sizeof(vertex_t)));
	out.area = nullptr;
	out.builtin_normals = 0;
	out.face_normals = nullptr;
	out.normals = nullptr;
	out.total_area = 0.0f;
	out.tree = nullptr;

	// Fill vertices
	vertex_t minPt; minPt.x = infOrMax<float>; minPt.y = infOrMax<float>; minPt.z = infOrMax<float>;
	vertex_t maxPt; maxPt.x = -infOrMax<float>; maxPt.y = -infOrMax<float>; maxPt.z = -infOrMax<float>;
	for (int i = 0; i < nverts; ++i) {
		vertex_t pt; pt.x = verts.at(i, 0, 0); pt.y = verts.at(i, 1, 0); pt.z = verts.at(i, 2, 0);
		minPt = vmin(minPt, pt);
		maxPt = vmax(maxPt, pt);
		out.vertices[i] = pt;
	}
	out.bBox[0] = minPt;
	out.bBox[1] = maxPt;

	// Fill faces
	for (int i = 0; i < nfaces; ++i) {
		face_t f; f.f0 = faces.at(i, 0, 0) - 1; f.f1 = faces.at(i, 1, 0) - 1; f.f2 = faces.at(i, 2, 0) - 1;
		out.faces[i] = f;
	}

	return out;
}

mxArray *makeFaceErrorStruct(const face_error *fe, int nfaces)
{
	static const char *faceErrorNames[7] = {
		"face_area",
		"min_error",
		"max_error",
		"mean_error",
		"mean_sqr_error",
		"serror",
		"samples_freq"
	};

	mxArray *mxFaceError = mxCreateStructMatrix(1, nfaces, 7, faceErrorNames);
	for (int i = 0; i < nfaces; ++i) {
		mxSetField(mxFaceError, i, "face_area", mxCreateDoubleScalar(fe[i].face_area));
		mxSetField(mxFaceError, i, "min_error", mxCreateDoubleScalar(fe[i].min_error));
		mxSetField(mxFaceError, i, "max_error", mxCreateDoubleScalar(fe[i].max_error));
		mxSetField(mxFaceError, i, "mean_error", mxCreateDoubleScalar(fe[i].mean_error));
		mxSetField(mxFaceError, i, "mean_sqr_error", mxCreateDoubleScalar(fe[i].mean_sqr_error));
		mxSetField(mxFaceError, i, "sample_freq", mxCreateDoubleScalar(fe[i].sample_freq));
		size_t nsamples = fe[i].sample_freq * (fe[i].sample_freq + 1) / 2;
		mxArray *mxSampleErrors = mxCreateDoubleMatrix(1, nsamples, mxREAL);
		double *sampleErrorData = mxGetPr(mxSampleErrors);
		for (int k = 0; k < nsamples; ++k) {
			sampleErrorData[k] = fe[i].serror[k];
		}
		mxSetField(mxFaceError, i, "serror", mxSampleErrors);
	}

	return mxFaceError;
}

mxArray *makeModelErrorStruct(const model_error& me)
{
	static const char *modelErrorNames[4] = {
		"min_error",
		"max_error",
		"mean_error",
		"n_samples"/*,
		"fe",
		"verror"*/
	};

	mxArray *mxModelError = mxCreateStructMatrix(1, 1, 4, modelErrorNames);
	mxSetField(mxModelError, 0, "min_error", mxCreateDoubleScalar(me.min_error));
	mxSetField(mxModelError, 0, "max_error", mxCreateDoubleScalar(me.max_error));
	mxSetField(mxModelError, 0, "mean_error", mxCreateDoubleScalar(me.mean_error));
	mxSetField(mxModelError, 0, "n_samples", mxCreateDoubleScalar(me.n_samples));
	/*mxSetField(mxModelError, 0, "fe", makeFaceErrorStruct(me.fe, me.mesh->num_faces));
	mxArray *mxVertError = mxCreateDoubleMatrix(1, me.mesh->num_vert, mxREAL);
	double *vertErrorData = mxGetPr(mxVertError);
	for (int i = 0; i < me.mesh->num_vert; ++i) {
		vertErrorData[i] = me.verror[i];
	}
	mxSetField(mxModelError, 0, "verror", mxVertError);*/
	return mxModelError;
}

mxArray *makeDistSurfSurfStatsStruct(const dist_surf_surf_stats& stats)
{
	static const char *distSurfSurfStatsNames[12] = {
		"st_m1_area",
		"m1_area",
		"m2_area",
		"min_dist",
		"max_dist",
		"mean_dist",
		"rms_dist",
		"cell_sz",
		"n_t_p_nec",
		"m1_samples",
		"grid_sz",
		"n_ne_cells"
	};

	mxArray *mxDistSurfSurfStats = mxCreateStructMatrix(1, 1, 12, distSurfSurfStatsNames);
	mxSetField(mxDistSurfSurfStats, 0, "st_m1_area", mxCreateDoubleScalar(stats.st_m1_area));
	mxSetField(mxDistSurfSurfStats, 0, "m1_area", mxCreateDoubleScalar(stats.m1_area));
	mxSetField(mxDistSurfSurfStats, 0, "m2_area", mxCreateDoubleScalar(stats.m2_area));
	mxSetField(mxDistSurfSurfStats, 0, "min_dist", mxCreateDoubleScalar(stats.min_dist));
	mxSetField(mxDistSurfSurfStats, 0, "max_dist", mxCreateDoubleScalar(stats.max_dist));
	mxSetField(mxDistSurfSurfStats, 0, "mean_dist", mxCreateDoubleScalar(stats.mean_dist));
	mxSetField(mxDistSurfSurfStats, 0, "rms_dist", mxCreateDoubleScalar(stats.rms_dist));
	mxSetField(mxDistSurfSurfStats, 0, "cell_sz", mxCreateDoubleScalar(stats.cell_sz));
	mxSetField(mxDistSurfSurfStats, 0, "n_t_p_nec", mxCreateDoubleScalar(stats.n_t_p_nec));
	mxSetField(mxDistSurfSurfStats, 0, "m1_samples", mxCreateDoubleScalar(stats.m1_samples));
	mxArray *mxGridSz = mxCreateDoubleMatrix(1, 3, mxREAL);
	double *gridSzData = mxGetPr(mxGridSz);
	gridSzData[0] = stats.grid_sz.x;
	gridSzData[1] = stats.grid_sz.y;
	gridSzData[2] = stats.grid_sz.z;
	mxSetField(mxDistSurfSurfStats, 0, "grid_sz", mxGridSz);
	mxSetField(mxDistSurfSurfStats, 0, "n_ne_cells", mxCreateDoubleScalar(stats.n_ne_cells));

	return mxDistSurfSurfStats;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//mexPrintf("Check inputs\n"); mexEvalString("drawnow;");
	ensureOrError(nrhs == 6, "Must supply 6 inputs");
	ensureOrError(isSize(prhs[0], { -1, 3 }), "FacesA must be an N x 3 array");
	ensureOrError(isSize(prhs[1], { -1, 3 }), "VerticesA must be an M x 3 array");
	ensureOrError(isSize(prhs[2], { -1, 3 }), "FacesB must be an N x 3 array");
	ensureOrError(isSize(prhs[3], { -1, 3 }), "VerticesB must be an M x 3 array");
	//mexPrintf("Get inputs\n"); mexEvalString("drawnow;");
	Volume<int> faceVolA = getVolumeChecked<int>(prhs[0], "FacesA");
	Volume<float> vertVolA = getVolumeChecked<float>(prhs[1], "VerticesA");
	Volume<int> faceVolB = getVolumeChecked<int>(prhs[2], "FacesB");
	Volume<float> vertVolB = getVolumeChecked<float>(prhs[3], "VerticesB");
	double sampleDensity = getCastScalarChecked<double>(prhs[4], "sampleDensity");
	int minSampleFreq = getCastScalarChecked<int>(prhs[5], "minSampleFreq");

	// Convert data to MESH data structures
	//mexPrintf("Init model errors\n"); mexEvalString("drawnow;");
	model_error me1 = initModelError();
	model_error me2 = initModelError();
	//mexPrintf("Build models\n"); mexEvalString("drawnow;");
	model m1 = buildModel(faceVolA, vertVolA);
	model m2 = buildModel(faceVolB, vertVolB);
	me1.mesh = &m1;
	me2.mesh = &m2;
	//mexPrintf("Compute distances\n"); mexEvalString("drawnow;");
	dist_surf_surf_stats stats1 = initDistSurfSurfStats();
	dist_surf_surf_stats stats2 = initDistSurfSurfStats();

	// Compute distance statistics
	//mexPrintf("Distance stats\n"); mexEvalString("drawnow;");
	dist_surf_surf(&me1, &m2, sampleDensity, minSampleFreq, &stats1, 0, nullptr);
	dist_surf_surf(&me2, &m1, sampleDensity, minSampleFreq, &stats2, 0, nullptr);

	int tmp1, tmp2;
	//mexPrintf("Vertex errors 1\n"); mexEvalString("drawnow;");
	calc_vertex_error(&me1, &tmp1, &tmp2);
	//mexPrintf("Vertex errors 2\n"); mexEvalString("drawnow;");
	calc_vertex_error(&me2, &tmp1, &tmp2);

	// Place results in MATLAB structs
	//mexPrintf("Put in plhs\n"); mexEvalString("drawnow;");
	plhs[0] = makeModelErrorStruct(me1);
	plhs[1] = makeModelErrorStruct(me2);
	/*plhs[2] = makeDistSurfSurfStatsStruct(stats1);
	plhs[3] = makeDistSurfSurfStatsStruct(stats2);*/

	// Free memory
	//mexPrintf("Free memory\n"); mexEvalString("drawnow;");
	free(me1.verror);
	free(me2.verror);
	free_face_error(me1.fe);
	free_face_error(me2.fe);
	mxFree(m1.faces);
	mxFree(m2.faces);
	mxFree(m1.vertices);
	mxFree(m2.vertices);

	//mexPrintf("Returning\n"); mexEvalString("drawnow;");
}
