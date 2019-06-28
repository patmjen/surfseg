% You may need to modify these to fit your system
GEL_LIB_ROOT_DIR = '../../../GEL';
MESH_LIB_ROOT_DIR = '../../../MESH';
BUILD_DIR = '../../build/surfseg/src';

COMMON_ARGS = {'matlab_util.cpp','-I../include','-I../src',...
    ['-I',GEL_LIB_ROOT_DIR,'/src'],...
    ['-I',MESH_LIB_ROOT_DIR],['-I',MESH_LIB_ROOT_DIR,'/lib3d/include'],...
    'CXXFLAGS="-std=c++14 -fPIC"',['-L',BUILD_DIR],'-lsurfseg'};

mex('mex_surfcut.cpp',COMMON_ARGS{:});
mex('mex_ksurfcut.cpp',COMMON_ARGS{:});
mex('mex_surfcut_planesep.cpp',COMMON_ARGS{:});
mex('mex_surfcut_planesep_dual.cpp',COMMON_ARGS{:});
mex('mex_surfcut_planesep_qpbo.cpp',COMMON_ARGS{:});

mex('mex_surfcut_4d.cpp',COMMON_ARGS{:});

mex('mex_gauss_curvature.cpp',COMMON_ARGS{:});
mex('mex_star_intersect.cpp',COMMON_ARGS{:});
mex('mex_subdiv_icosahedron.cpp',COMMON_ARGS{:});
mex('mex_hausdorff.cpp',COMMON_ARGS{:});
