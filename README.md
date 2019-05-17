# gds2Para
Complete Integrated Circuit (IC) Layout Analysis from GDSII Design File to Parasitics Extraction

This layout analyzer is written in C++ as part of a wider API for the electromagnetic design validation of VLSI IC, package, and board designs. The code design follows a philosophy intended to have the most extensibility and highest level of automation possible for increasingly complex integrated circuits, packages, and boards. Design information from the Graphic Database Stream II (GDSII) file is parsed and stored alongside simulation input information. Once a design is loaded, the information can be written to another GDSII file or analyzed with a full-wave simulation to extract parasitics between ports. An alternate mode of operation allows an interconnect modeling platform (IMP) file to be translated to a GDSII file design. Parasitics are preferably reported through a [Xyce](https://xyce.sandia.gov/) subcircuit file (SPICE derivative from Sandia National Laboratories), though a [Standard Parasitic Exchange Format (SPEF)](https://ieeexplore.ieee.org/document/5430852) file may be output alternatively. Please carefully read this **README.md** file for an overview, dependencies, and usage instructions.

## Overview
### Simulate (-s) Mode Top-level Flowchart
<img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/gds2Para_mode-s.png" width=600 alt="Simulate (-s) Mode Flowchart">

## Packages
| Packages / Library               | Languages        | Headers                                | Description                                                       |
| -------------------------------- | ---------------- | -------------------------------------- | ----------------------------------------------------------------- |
| C Standard Library               | C                | CCtype, Cerrno, Cmath, Cstdio, Cstlib, Cstring, Ctime     | Utilities for testing character categorization, C language error numbers, basic arithmetic operations, C language input/output, manipulating the char\* data type, and dates and timers |
| C++ Standard Library: Algorithms | C++              | Algorithm                              | Useful algorithms including sorting                               |
| C++ Standard Library: Containers | C++              | Queue, Set, Stack, Unordered\_map, Unordered\_set, Vector | Utilities for handling queues and priority queues, sets and multisets, LIFO stacks, associative data types, dictionary collection types, and collection data types from existing C++ data types  |
| C++ Standard Library: I/O        | C++              | FStream, IOStream                      | Plain text file input/output and command line input/output        |
| C++ Standard Library: Numerics   | C++              | Complex                                | Complex floating point number arithmetic                          |
| C++ Standard Library: Strings    | C++              | String                                 | Utilities for easy string handling                                |
| C++ Standard Library: Utilities  | C++              | Utility                                | General C++ utilities such as pairs and swapping                  |
| [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)                                                             | C++              | Sparse                                 | API for sparse matrix storage                                     |
| [HYPRE](https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods)                            | Fortran/C/C++    | HYPRE, HYPRE_krylov, HYPRE_parcsr_ls   | LLNL library of preconditioners and multigrid solvers             |
| Limbo                            | C/C++            | GdsReader, GdsWriter                   | API for GDSII file parsing and writing                            |
| [MakeUtils](https://www.gnu.org/software/make/)                                                                           | Makefile         | n/a                                    | Makefile utilities that help find build dependencies              |
| [Math Kernel Library](https://software.intel.com/en-us/mkl)                                                               | Fortran/C/C++    | MKL, Mkl\_spblas                        | Intel libraries for parallel mathematics                          |
| [Open Message Passing Interface](https://www.open-mpi.org/)                                                               | Fortran/C        | n/a                                    | Open MPI compilation and runtime tools for parallel computing     |
| Parser-SPEF                      | C++              | parser-spef                            | API for SPEF file processing                                      |
| Preprocessor                     | C/C++            | n/a                                    | Inclusion, definitions, and macros such as assertion with message |

Note that modern compilers that support the following language standards are needed to make the package dependencies: GNU Fortran (superset of F95) for Fortran, C99 for C language, and C++17 for C++.

## Custom Classes
| Structure / Class Name   | Header                  | Description                                                       |
| ------------------------ | ----------------------- | ------------------------------------------------------------------|
| AsciiDataBase            | limboint                | Full GDSII design being read/constructed with file name, metadata, units, data of all geometric cells, and conductor information |
| GeoCell                  | limboint                | GDSII geometric cell with name, metadata, and all GDSII elements  |
| boundary                 | limboint                | GDSII element for a polygonal boundary                            |
| path                     | limboint                | GDSII element for a path                                          |
| node                     | limboint                | GDSII element for an electrical node                              |
| box                      | limboint                | GDSII element for a box outline                                   |
| textbox                  | limboint                | GDSII element for a text box                                      |
| sref                     | limboint                | GDSII element for a structure reference                           |
| SolverDataBase           | solnoutclass            | Design name, output files, information other than layout, and all results of the simulation necessary for parameter extraction   |
| SimSettings              | solnoutclass            | Simulation settings including units, design limits, and frequency sweep parameters                                               |
| Layer                    | solnoutclass            | Layer information from physical stack-up including name, GDSII layer number, z-coordinates, and material properties              |
| Waveforms                | solnoutclass            | Needed information for time-domain plots (unimplemented)          |
| Parasitics               | solnoutclass            | Collection of ports and sparse matrices for admittance parameters |
| Port                     | solnoutclass            | Port information including name, reference direction, supply and return coordinates, and containing GDSII layer                  |
| fdtdMesh                 | fdtd                    | Master class for all discretization, conductor, and solver data   |
| fdtdPatch                | fdtd                    | Information of a single discretized patch                         |
| fdtdBound                | fdtd                    | Information of a single conductor or dielectric boundary region   |
| fdtdCdt                  | fdtd                    | Solver information common to all conductors                       |
| fdtdOneCondct            | fdtd                    | Solver information regarding discretization and excitation of a single conductor region                                          |
| fdtdPort                 | fdtd                    | Solver representation of a port                                   |

## Developers
All software is under active development without any notices regarding the timing of nature of updates. Contributions are welcome through forking and [pull requests](https://github.com/purdue-onchip/gds2Para/pulls). Please direct all questions, bug reports, and issues relating to this repository to the primary maintainer:
* [Michael Hayashi](mailto:mhayashi@purdue.edu?subject=Inquiry%20for%20gdsii-interface), Graduate Research Assistant, School of ECE, Purdue University

Developers:
* Michael Hayashi, Graduate Research Assistant, School of ECE, Purdue University
* Li Xue, Graduate Research Assistant, School of ECE, Purdue University

Many thanks are owed to the project sponsors for making the development of this software possible:

| <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/purdue.png" width=100 alt="Purdue University"> | <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/DARPA.png" width=100 alt="DARPA ERI"> |
| :---: | :---: |

## Installation and Usage
Follow the instructions given in **INSTALL.md** to install this software and run it from the command line. Users from Purdue University should read **purdue\_install.md** instead for specific steps unique to their environment.

## Credits and Acknowledgements
* Example Files
    * [Purdue On-Chip Electromagnetics Laboratory](https://engineering.purdue.edu/~djiao/) for **singleStrip.imp** (single stripline)
    * [Yuan Ze University Design Automation library](http://yzuda.org/download/_GDSII_examples.html) for **inv.gds** (inverter logic gate), **nand2.gds** (2-input NAND gate), and **xor.gds** (XOR gate)
    * [NanGate Open Cell Library FreePDK45](https://www.cs.upc.edu/%7Ejpetit/CellRouting/nangate/Front_End/Doc/Databook/Cells/SDFFRS_X2_NangateOpenCellLibrary_typical_typical.html) for **SDFFRS\_X2.gds** (scan D flip-flop with reset/set)
    * [Peter Monta](https://github.com/pmonta/FPGA-netlist-tools/tree/master/masks/4004) and the [Tuva Design Intel 4004 Anniversary Project](http://4004.com/) for **4004.gds** (Intel 4004 processor from 1971)
    * [Open Compute Project](https://www.opencompute.org/wiki/Networking/SpecsAndDesigns#Facebook_Wedge_100_-_32x100G) for original Gerber files for Facebook Wedge 100 that led to **WEDGE100_REV01.gds** (100G Ethernet Switch board)
* Prof. Dan Jiao of the Purdue On-Chip Electromagnetics Library for her guidance of the project
* The staff of the Defense Advanced Research Projects Agency (DARPA) Microsystems Technology Office (MTO) for their sponsorship, feedback, and support of this endeavor
* The other contributors of the [Intelligent Design of Electronic Assets (IDEA)](https://github.com/aolofsson/IDEA) team and the [Posh Open Source Hardware (POSH)](https://github.com/aolofsson/POSH) team within the DARPA Electronics Resurgence Initiative (ERI)

## License
GNU Public License v2.0 (see **LICENSE**)
