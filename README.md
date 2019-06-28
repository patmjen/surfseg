# SurfSeg
Volumetric segmentation using explicit surfaces.

This is a small C++14 library to perform segmentation in volumetric images with explicit (hyper)surfaces defined via meshes.
A C++ API is provided as well as a set of mex files for MATLAB.

Some MATLAB scripts have been created to show how to use the library.

## How to build
The project can be built using [CMake](https://cmake.org/). To build, navigate to the project directory and run:
```
mkdir build
cd build
cmake ..
```
and then run the resulting generator (most likely a makefile or Visual Studio project).

In order to build this project you will need the [GEL](https://github.com/janba/GEL) library.
If you use Windows and compile GEL with provided Visual Studio files, you can set the environment variable `GEL_LIB_ROOT_DIR`
before calling CMake to specify where GEL is stored. However, you may have to tweak this manually.

You will also need the [MESH](https://github.com/naspert/MESH) tool to compile the mex function which computed the Hausdorff
distance between two meshes. Again, you can set the `MESH_LIB_ROOT_DIR` environment variable before calling CMake to specify
where the tool is stored -- but you will likely need to adjust the CMakeLists.txt manually to get linking to work.
**If you want to skip this** and don't care about this mex function, then set the CMake variable `SKIP_MEX_HAUSDORFF` to true.