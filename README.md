# gds2Para
Complete Integrated Circuit (IC) Layout Analysis from GDSII Design File to Parasitics Extraction

This IC layout analyzer is written in C++ as part of a wider API for the electromagnetic design validation of VLSI designs. The code design follows a philosophy intended to have the most extensibility and highest level of automation possible for increasingly complex integrated circuits. Design information from the Graphic Database Stream II (GDSII) file is parsed and stored alongside simulation input information. Once a design is loaded, the information can be written to another GDSII file or analyzed with a full-wave simulation to extract parasitics between ports. An alternate mode of operation allows an interconnect modeling platform (IMP) file to be translated to a GDSII file design. Parasitics are reported through a Standard Parasitics Exchange Format (SPEF) file. Please carefully read this **README.md** file for installation and usage instructions.

## Overview
### Simulate (-s) Mode Top-level Flowchart
| <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/gds2Para_mode-s.png" width=100 alt="Simulate (-s) Mode Flowchart">

## Packages
| Packages                                | Languages                       | Description                                                                |
| --------------------------------------- | ------------------------------- | -------------------------------------------------------------------------- |
| Algorithms                              | C++                             | Useful algorithms including sorting                   |
| Cctype                                  | C                               | Utilities for testing character categorization     |
| Cerrno                                  | C                               | Utilities for C language error numbers             |
| Cmath                                   | C                               | Utilities for basic arithmetic operations          |
| Complex                                 | C++                             | Complex number arithmetic                         |
| Cstdio                                  | C                               | Utilities for C language input/output              |
| Cstdlib                                 | C                               | Utilities for common C language tasks              |
| Cstring                                 | C                               | Utilities for manipulating char\* data type         |
| Ctime                                   | C                               | Utilities for dates and timers                      |
| Eigen                                   | C++                             | API for sparse matrix storage                       |
| FStream                                 | C++                             | Plain textext file input and output                        |
| GDS2GDT                                 | C                               | Binary-to-text conversion from .gds to .gdt file           |
| GeoCell                                 | C++                             | Custom class for GDSII geometric cell             |
| Geometry                                | C++                             | Custom class for VLSI design geometry             |
| IOStream                                | C++                             | Command line input and output                     |
| Limbo                                   | C/C++                           | API for GDSII file parsing and writing                      |
| MakeUtils                               | Makefile                        | Makefile utilities that help find dependencies                             |
| Math                                    | C++                             | Extension of math functions from STL                                       |
| MKL                                     | Fortran/C/C++                   | Intel libraries for parallel mathematics          |
| Parser-SPEF                             | C++                             | API for SPEF file processing                      |
| Preprocessor                            | C/C++                           | Macros such as assertion with message                                      |
| ProgramOptions                          | C++                             | Easy API to parser command line options                                    |
| Set                                     | C++                             | Utilities for sets and multisets                  |
| Stack                                   | C++                             | Utilities for handling LIFO stacks                 |
| String                                  | C++                             | Utilities for easy string handling                                              |
| ThirdParty                              | C/C++                           | Third party packages required                                              |
| Unordered\_map                           | C++                             | Utilities for handling associative data types     |
| Unordered\_set                           | C++                             | Utilities for handling dictionary collection types   |
| Utility                                 | C++                             | Functions and classes for pairs and swapping       |
| Vector                                  | C++                             | Collection data types from existing C++ data types |

## Developers
Please direct all questions, bug reports, and issues relating to this repository to the primary maintainer:
* [Michael Hayashi](mailto:mhayashi@purdue.edu?subject=Inquiry%20for%20gdsii-interface), Graduate Research Assistant, School of ECE, Purdue University

Developers:
* Michael Hayashi, Graduate Research Assistant, School of ECE, Purdue University
* Li Xue, Graduate Research Assistant School of ECE, Purdue University

Many thanks are owed to the project sponsors:

| <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/purdue.png" width=100 alt="Purdue University"> | <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/DARPA.png" width=100 alt="DARPA ERI"> |
| :---: | :---: |

## Installation and Usage
1. Clone this respository into your working location
2. Repeat for the [Limbo repository](https://github.com/limbo018/Limbo), [Parser-SPEF repository](https://github.com/OpenTimer/Parser-SPEF), and [eigen-git-mirror repository](https://github.com/eigenteam/eigen-git-mirror)
3. Ensure that you are the owner of the files that were downloaded with full read/write/execute permissions
4. Ensure that GNU bison is installed on the machine by seeing if there is an output to the terminal with `which bison`
5. Modify run commands files depending on the shell indicated by `echo $SHELL`, substituting \<absolute path to Limbo directory> for a valid path when it appears
6. Exit the shell and terminate the connection before logging back in
7. Enter the directory "Limbo/limbo/parsers/gdsii/stream" from the working location
8. Run `make CXX=g++ CC=gcc FC=gfortran` in this directory
9. Return to working directory and ensure that a new library file named **libgdsparser.a** exists by running `ls -lh Limbo/lib`
10. Enter the directory "gds2Para" from the working location
11. (Optional) Add GDSII file to working location
12. Run `make` in shell to compile executable **LayoutAnalyzer**
13. Run `LayoutAnalyzer --help` in shell to view possible modes of operation for this software

## License
GNU Public License v2.0 (see **LICENSE**)
