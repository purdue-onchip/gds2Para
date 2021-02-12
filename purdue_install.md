# Instructions Specific to Purdue ECN Machines
Please follow these instructions for a more complete deployment on Red Hat Enterprise Linux (RHEL) machines at Purdue University maintained by the Engineering Computer Network (ECN). There are a number of nonuniformities that can cause subtle, pernicious failures when attempting to deploy this code for yourself, but this guide attempts to address them as they are discovered. These instructions do **not** apply to any users outside Purdue University.

## Developers
Please direct all questions, bug reports, and issues relating to installing or running this software at Purdue University to the primary maintainer:
* [Michael Hayashi](mailto:mhayashi@purdue.edu?subject=Inquiry%20for%20gds2Para), Graduate Research Assistant, School of ECE, Purdue University

## Installation and Usage
1. Ensure that the following programs are installed on the machine by seeing if there is an output to each of the following commands in these categories:
    * Selected package manager development tools
        * GNU Autoconf: `which autoconf`
        * GNU Automake: `which automake`
        * GNU Bison: `which bison`
        * Flex: `which flex`
        * GNU Compiler Collection (GCC version 7 or newer): `which gcc`, `which g++`, and `which gfortran`
        * GNU Make: `which make`
    * Boost C++ Libraries (Boost version 1.66.0 or newer): `locate "filtering_stream.hpp"`
    * CMake (CMake version 3.7.0 or newer): `which cmake` and `cmake --version`
    * git: `which git`
    <!--* Intel Math Kernel Library (MKL): `locate "mklvars.sh"`
    * Lemon: `which lemon`
    * OpenMPI (OpenMPI version 4.0.1 or newer): `which mpicc`, `which mpicxx` -->
    * Zlib (Zlib version 3.13 or newer): `locate "zlib.h"`
    * If any of these programs are not installed, then it is recommended that someone with root access install whatever is missing using the package manager (most programs) or by building from source (GCC, Boost, and OpenMPI).
2. Clone this respository into your working location with `git clone git@github.com:purdue-onchip/gds2Para.git`
    * If cloning fails, the usual reason is that a secure connection to the server could not be made. It is important to clone to be able to use `git` to receive updates. These steps will help get yourself going using `git`:
        1. Check if you have any public keys such as **id_rsa.pub** by running `ls -lh .ssh` in your home directory
        2. If no public keys appear, then generate one by running `ssh-keygen -t rsa -C <Purdue email>` in your home directory, where \<Purdue email> is your Purdue email address ending in "@purdue.edu"
        3. The key generation command will require you to click "Enter" to place the public key in the default location, and then click "Enter" twice more to forgo additional password protection
        4. Open the public key **id_rsa.pub** in a text editor; for example, run `vim ~/.ssh/id_rsa.pub`
        5. In your GitHub account, go to [your SSH and GPG keys under your account settings](https://github.com/settings/keys) and click the green "New SSH key" button
        6. Give this SSH public key title such as "Purdue ECN"
        7. Copy and paste the contents of **id_rsa.pub** open in the text editor verbatim into the key field
        8. Click the green "Add SSH key" button
        9. Navigate back to your working location and attempt to clone the repository again
    * If cloning and SSH key generation fail, try to download the files by clicking the green "Clone or download" button on the Code tab of this repository followed by "Download ZIP"
        1. Move the zip archive to the working location
        2. Unzip the archive with `unzip gds2Para-master.zip`
        3. Rename the recently created directory from "gds2Para-master" to "gds2Para"
3. Clone the [Limbo repository](https://github.com/limbo018/Limbo) into your working location with `git clone https://github.com/limbo018/Limbo.git`
    * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Limbo repository followed by "Download ZIP"
        1. Move the zip archive to the working location
        2. Unzip the archive with `unzip Limbo-master.zip`
        3. Rename the recently created directory from "Limbo-master" to "Limbo"
4. Clone the [Parser-SPEF repository](https://github.com/OpenTimer/Parser-SPEF) into your working location with `git clone https://github.com/OpenTimer/Parser-SPEF.git`
    * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Parser-SPEF repository followed by "Download ZIP"
        1. Move the zip archive to the working location
        2. Unzip the archive with `unzip Parser-SPEF-master.zip`
        3. Rename the recently created directory from "Parser-SPEF-master" to "Parser-SPEF"
5. Clone the [eigen-git-mirror repository](https://github.com/eigenteam/eigen-git-mirror) into your working location with `git clone https://github.com/eigenteam/eigen-git-mirror.git`
    * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the eigen-git-mirror repository followed by "Download ZIP"
        1. Move the zip archive to the working location
        2. Unzip the archive with `unzip eigen-git-mirror-master.zip`
        3. Rename the recently created directory from "eigen-git-mirror-master" to "eigen-git-mirror"
6. Ensure that you are the owner of the files that were downloaded with full read/write/execute permissions
    * Change ownership with `chown -R <username>:<username> <directory>`, where \<username> is your Purdue username and \<directory> is each of "gds2Para/", "Limbo/", "Parser-SPEF/", or "eigen-git-mirror/"
    * Change permissions with `chmod -R 744 <directory>`, where \<directory> is each of "gds2Para/", "Limbo/", "Parser-SPEF/", or "eigen-git-mirror/"
7. Modify run commands files depending on the shell indicated by `echo $SHELL`, substituting \<absolute path to Limbo directory>, \<absolute path to Parser-SPEF directory>, \<absolute path to Eigen directory>, \<absolute path to future HYPRE directory> for valid paths when they appears
    * For "bash", edit the file ".bashrc" in your home directory by appending the following:
    ```bash
    # Skip rest of file if not interactive
    if [ -z "$PS1" ]; then
        return
    fi

    # Note: Commands for improving "gcc" and "OpenMPI" can be replaced by 
    #       command "module load gcc" and "module load mpi" in new OS.
    # Details in https://engineering.purdue.edu/ECN/Support/KB/Docs/LmodonLinuxHosts
    module load gcc/8.3.0
    module load mpi/mpich-3.2-x86_64

    # Improve Compilers
    alias cc="gcc"
    export CC=/package/gcc/8.3.0/bin/gcc
    export CXX=/package/gcc/8.3.0/bin/g++
    export FC=/package/gcc/8.3.0/bin/gfortran
    export PATH=/package/gcc/8.3.0/bin:$PATH
    export LD_LIBRARY_PATH=/package/gcc/8.3.0/lib64:/package/gcc/8.3.0/lib:$LD_LIBRARY_PATH

    # Improve OpenMPI Version
    export PATH=/usr/lib64/mpich-3.2/bin:${PATH}
    export LD_LIBRARY_PATH=/usr/lib64/mpich-3.2/lib:${LD_LIBRARY_PATH}
    export OMPI_CC=$CC
    export OMPI_CXX=$CXX
    export OMPI_FC=$FC

    # Improve Intel MKL Version (module load intel)
    module load intel/18.0.1
    export MKL_DIR=/package/intel/18.0.1/mkl
    export MKLROOT=$MKL_DIR
    export LD_LIBRARY_PATH=$MKL_DIR/lib/intel64:$LD_LIBRARY_PATH

    # Improve git Version
    export PATH=/opt/git/2.18.0/bin:$PATH
    export LD_LIBRARY_PATH=/opt/git/2.18.0/lib64:/opt/git/2.18.0/libexec:$LD_LIBRARY_PATH

    # Environment Variables for GDSII File Handling
    export BOOST_DIR="/opt/boost/1.57.0"
    export LIMBO_DIR="<absolute path to Limbo directory>"
    export PARSER_SPEF_DIR="<absolute path to Parser-SPEF directory>"
    export EIGEN_DIR="<absolute path to Eigen directory>"
    export HYPRE_DIR="<absolute path to future HYPRE directory>"
    ```
    * For "tcsh", edit the file ".cshrc" in your home directory by appending the following:
    ```tcsh
    # Skip Rest of File if Not Interactive
    if (! $?prompt) then
        exit 0
    endif

    # Note: Commands for improving "gcc" and "OpenMPI" can be replaced by 
    #       command "module load gcc" and "module load mpi" in new OS.
    # Details in https://engineering.purdue.edu/ECN/Support/KB/Docs/LmodonLinuxHosts
    module load gcc/8.3.0
    module load mpi/mpich-3.2-x86_64

    # Improve gcc Version
    alias cc 'gcc'
    setenv CC /package/gcc/8.3.0/bin/gcc
    setenv CXX /package/gcc/8.3.0/bin/g++
    setenv FC /package/gcc/8.3.0/bin/gfortran
    setenv PATH /package/gcc/8.3.0/bin:$PATH
    setenv LD_LIBRARY_PATH /package/gcc/8.3.0/lib64:/package/gcc/8.3.0/lib:${LD_LIBRARY_PATH}

    # Improve OpenMPI Version
    setenv PATH /usr/lib64/mpich-3.2/bin:${PATH}
    setenv LD_LIBRARY_PATH /usr/lib64/mpich-3.2/lib:${LD_LIBRARY_PATH}
    setenv OMPI_CC ${CC}
    setenv OMPI_CXX ${CXX}
    setenv OMPI_FC ${FC}

    # Improve Intel MKL Version
    module load intel/18.0.1
    setenv MKL_DIR /package/intel/18.0.1/mkl
    setenv MKLROOT ${MKL_DIR}
    setenv LD_LIBRARY_PATH ${MKL_DIR}/lib/intel64:${LD_LIBRARY_PATH}

    # Improve git Version
    setenv PATH /opt/git/2.18.0/bin:${PATH}
    setenv LD_LIBRARY_PATH /opt/git/2.18.0/lib64:/opt/git/2.18.0/libexec:${LD_LIBRARY_PATH}

    # Environment Variables for GDSII File Handling
    setenv BOOST_DIR /opt/boost/1.57.0
    setenv LIMBO_DIR <absolute path to Limbo directory>
    setenv PARSER_SPEF_DIR <absolute path to Parser-SPEF directory>
    setenv EIGEN_DIR <absolute path to Eigen directory>
    setenv HYPRE_DIR <absolute path to future HYPRE directory>
    ```
8. Exit the shell and terminate the connection before logging back in
9. Ensure that the run command files were properly loaded by running `echo $LIMBO_DIR`
    * If nothing shows up for "bash" users, run `cp .bashrc .bash_profile` in the home directory, exit the shell, log back in, and try again
    * If nothing shows up for "tcsh" users, run `cp .cshrc .tcshrc` in the home directory, exit the shell, log back in, and try again
    * For all other errors, contact the primary maintainer
10. Enter the directory "Limbo" from the working location
11. Run `git reset --hard 3.3.0` to revert to an earlier version of the repository
    * The partial installation of Limbo components needed to run gds2Para is feasible only before CMake was introduced as the way to build the entirety of Limbo
    * Only the necessary parsers will be built from Limbo using these instructions, and users wanting to use other Limbo features should follow the installation instructions given in that repository
    * (Optional) It is believed that the necessary components and most other components of Limbo can be built using `make` in the top directory of Limbo following these instructions, but this has not been tested
12. Install third-party LEF & DEF readers by following the instructions in file /Limbo/limbo/thirdparty/lefdef/5.8/lefdefReadme.txt
13. Compile Limbo's parsers 
    * Limbo's GDSII parser: Enter the directory "/Limbo/limbo/parsers/gdsii/stream". Run `make -j 2` in this directory to build certain Limbo libraries
    * Limbo's LEF parser: Enter the directory "/Limbo/limbo/parsers/lef/adapt". Run `make -j 2`.
    * Limbo's DEF parser: Enter the directory "/Limbo/limbo/parsers/def/adapt". Run `make -j 2`.
14. Ensure that new static libraries named **libgdsparser.a**, **liblefparseradapt.a**, and **libdefparseradapt.a** exist in directory "/Limbo/lib"
15. Follow the instructions in [HYPRE Setup](#HYPRE-Setup) to prepare multigrid methods from HYPRE
    * It is necessary to make all changes to ".bashrc" and ".cshrc" in the home directory before following these instructions
    * Continue following the remaining steps once HYPRE is installed
16. Enter the directory "gds2Para" from the working location
17. (Optional) Add GDSII files besides "nand2.gds", "SDFFRS_X2.gds", and "4004.gds" to the directory "examples/" within the working location
    * Some usage modes of this software may also require a simulation input file in a similar format to "nand2.sim_input", "SDFFRS_X2.sim_input", and "4004.sim_input" (the file name before the suffix must match that of the GDSII file before ".gds")
    * Files created from the software will be added directly to the present working directory, which is "gds2Para" for this point in these instructions
18. Run `make` in shell to compile executable **LayoutAnalyzer**
19. Run `LayoutAnalyzer --help` to get a list of the available control modes which the executable supports
20. Run `LayoutAnalyzer -r examples/nand2.gds` in shell to produce terminal output describing the GDSII file ("nand2.gds" may be replaced with any GDSII file available)
21. Perform a complete parameter extraction by running `mpirun LayoutAnalyzer -s examples/SDFFRS_X2.gds examples/SDFFRS_X2.sim_input examples/SDFFRS_X2.cir` to read in the design and simulation input file, do all analysis, and return the results in a Xyce (SPICE-compatible) subcircuit
    * It is necessary to use `mpirun` to guarantee memory integrity in the parallelized portions of the software.

## HYPRE Setup
1. Clone the [HYPRE repository](https://github.com/hypre-space/hypre) into your working location with `git clone https://github.com/hypre-space/hypre.git`
    * This requires the setup of an SSH key for a secure connection
    * Setup will follow the instructions at "hyper/INSTALL" using information known to affect Purdue ECN users
2. Enter the directory "hypre/src/cmbuild" from the working location
    * Configuration options may be viewed using `../configure --help`
<!--3. Run the command `./configure --enable-shared --enable-complex --with-MPI=/opt/mpich2-gcc/1.4.1p1 --with-MPI-include=/opt/mpich2-gcc/1.4.1p1/include/ --with-MPI-libs="" --with-MPI-lib-dirs=/opt/mpich2-gcc/1.4.1p1/lib64/ CC=/opt/gcc/7.1.0/bin/gcc CFLAGS="-g" CXX=/opt/gcc/7.1.0/bin/gcc CXXFLAGS="-g" FC=/opt/gcc/7.1.0/bin/gfortran`-->
3. Run `cmake --version` and ensure that the version of CMake installed is at least 2.8.8
4. Run `cmake -LA ..` and ensure that CMake found the correct versions of the compilers
    * CMake is rather temperamental and tends to overlook common methods of specifying different compilers system-wide
    * Check ".bashrc" or ".cshrc" in your home directory to see if the environment variables for the compilers are properly set
    * For all other issues, contact the primary maintainer
5. Run `cmake ..` to produce a makefile
6. Run `make install -j 2` from within the same location ("hypre/src/cmbuild" in the working location) to install header files
7. Return to working directory and ensure that a new library file named **HYPRE.h** exists by running `ls -lh hypre/src/hypre/include`
