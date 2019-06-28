mex mex_surfcut.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_ksurfcut.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_surfcut_planesep.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_surfcut_planesep_dual.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_surfcut_planesep_qpbo.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg

mex mex_surfcut_4d.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg

mex mex_gauss_curvature.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_star_intersect.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_subdiv_icosahedron.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src CXXFLAGS="-std=c++14 -fPIC" -L../../build/surfseg/src -lsurfseg
mex mex_hausdorff.cpp matlab_util.cpp -I../include -I../src -I../../../GEL/src -I../../../MESH -I../../../MESH/lib3d/include CXXFLAGS="-std=c++14 -fPIC" -LDFLAGS="-fPIC" -L../../build/surfseg/src -L../../../MESH/build -lsurfseg -lmesh
