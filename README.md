# gdsii-interface
GDSII File Parser and Plotter

This GDSII file parser is written in C++ as part of a wider API for the electromagnetic design validation of VLSI designs. The code design follows a philosophy intended to have the most extensibility and highest level of automation possible for increasingly complex integrated circuits. Please carefully read this **README.md** file for installation and usage instructions

## Packages
| Packages                                | Languages                       | Description                                                                |
| --------------------------------------- | ------------------------------- | -------------------------------------------------------------------------- |
| Algorithms                              | C++                             | Useful algorithms including partitioning, coloring, etc.                   |
| FStream                                 | C++                             | Text file input and output                        |
| GDS2GDT                                 | C                               | Binary-to-text conversion from .gds to .gdt file           |
| GeoCell                                 | C++                             | Custom class for GDSII geometric cell             |
| Geometry                                | C++                             | Custom class for VLSI design geometry             |
| gnuplot                                 | Shell                           | Command-line utility for plotting                 |
| IOStream                                | C++                             | Command line input and output                     |
| MakeUtils                               | Makefile                        | Makefile utilities that help find dependencies                             |
| Math                                    | C++                             | Extension of math functions from STL                                       |
| Preprocessor                            | C++                             | Macros such as assertion with message                                      |
| ProgramOptions                          | C++                             | Easy API to parser command line options                                    |
| String                                  | C++                             | Utilities to char and string                                               |
| ThirdParty                              | C/C++                           | Third party packages required                                              |
| Vector                                  | C++                             | Collection data types from existing C++ data types |

## Developers
Please direct all questions, bug reports, and issues relating to this repository to the primary maintainer:
* [Michael Hayashi](mailto:mhayashi@purdue.edu?subject=Inquiry%20for%20gdsii-interface), Graduate Research Assistant, School of ECE, Purdue University

Many thanks are owed to the project sponsors:

| <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/purdue.png" width=100 alt="Purdue University"> | <img src="https://github.com/purdue-onchip/gdsii-interface/blob/master/images/DARPA.png" width=100 alt="DARPA ERI"> |
| :---: | :---: |

## Installation and Usage
1. Clone this respository into your working location
2. Enter "gdsii-interface" directory in working location
3. (Optional) Add GDSII file to working location
4. Run `make all` in shell to compile executable **TestInterface**
5. (Optional) Modify **gdsplot.sh** to work on GDSII file added
6. Run `bash ./gdsplot.sh` in shell to produce PNG images recreating each layer in the GDSII file

## License
GNU Public License v2.0 (see **LICENSE**)
