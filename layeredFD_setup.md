# How to Setup Layered Finite-Difference Solver
Please follow these instructions for running layered FD solver in Linux system or Windows system. Please direct all questions, bug reports, and issues to the primary maintainer:
* [Shuzhan Sun](mailto:sun630@purdue.edu?subject=Inquiry%20for%20gds2Para), Graduate Research Assistant, School of ECE, Purdue University

## Linux System
The macro `SKIP_LAYERED_FD` in head file `fdtd.hpp` switches between layeredFD solver and V0Vh solver. Setup for each is:
* Case 1: run V0Vh solver in Linux
```C++
#define SKIP_LAYERED_FD
```

* Case 2: run layeredFD solver in Linux
```C++
// #define SKIP_LAYERED_FD
```

Everything else is the same as [V0Vh solver](https://github.com/purdue-onchip/gds2Para/blob/master/purdue_install.md), including commands, final parameters, and storage styles.

## Windows System
Currently, only layeredFD solver is supported in Windows system, and above macro `SKIP_LAYERED_FD` no longer changes anything in Windows system.
### How to Setup in Visual Studio (tested in Visual Studio 2017):
1. Create a visual studio project next to folder `gds2Para`, and only include these files in project:
```bash
fdtd.hpp
generateStiff.cpp
layeredFdtd.hpp
mapIndex.hpp
matrixTypeDef.hpp
mesh.cpp
pardisoSolver.hpp
sysInfoIO.hpp
```
2. Edit the `configuration` of the visual studio project to include Intel MKL. A simple way is to install [Intel Parallel Studio XE](https://software.intel.com/en-us/parallel-studio-xe/choose-download), which is free for student.
3. Write a main function to call the layeredFD solver and put it next to the VS project solution (.sln). Example main function is `layeredFD.cpp`:
```C++
#include "gds2Para/src/layeredFdtd.hpp"
int main(void) {
    layeredFdtd();
    return 0;
}
```
### How to Run in Visual Studio
1. Run the layeredFD solver in Linux (could stop ealier with `Ctrl+Z`) to generate a folder `temp_sysInfoIO` next to folder `gds2Para`. This step exports the structure information to a few txt files. This step is necessary because many codes for loading gds files cannot be compiled in Windows system.
2. Back to Visual Studio, under mode `Debug/Release x64/x86`, run `Local Windows Debugger` in VS to get the Z-parameters.
