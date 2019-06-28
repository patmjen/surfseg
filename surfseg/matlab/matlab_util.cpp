#include "matlab_util.h"
#include <algorithm>

bool isSize(const mxArray *a, std::initializer_list<int> size)
{
	size_t dimA = mxGetNumberOfDimensions(a);
	int i = 0;
	for (int s : size) {
		if (s != -1 && (i >= dimA || s != mxGetDimensions(a)[i])) {
			return false;
		}
		++i;
	}
	return true;
}

int getMaxCompThreads()
{
	mxArray *matlabCallOut[1] = { 0 };
	mxArray *matlabCallIn[1] = { 0 };
	mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
	double *Nthreadsd = mxGetPr(matlabCallOut[0]);
	return (int)Nthreadsd[0];
}
