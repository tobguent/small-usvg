# small-usvg #

This folder contains a reference implementation of the paper:

> *Xingze Tian, Tobias Günther*  
**Unified Smooth Vector Graphics: Modeling Gradient Meshes and Curve-based Approaches Jointly as Poisson Problem**  
IEEE Transactions on Visualization and Computer Graphics

The implementation was tested on MSVC and GCC. 
A CMake file is provided to compile the program.

Folders:
- `usvg` *Contains a C++ library that can read scenes, build patches, and solve PDEs.*
- `demo` *Calls the library to compute the result images of the paper, which are written into the build folder.*