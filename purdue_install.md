# Instructions Specific to Purdue ECN Machines
Please follow these instructions for a more complete deployment on Red Hat Enterprise Linux (RHEL) machines at Purdue University maintained by the Engineering Computer Network (ECN). There are a number of nonuniformities that cause subtle, pernicious failures when attempting to deploy this code for yourself, but this guide attempts to address them as they are discovered. These instructions do **not** apply to any users outside Purdue University.

## Developers
Please direct all questions, bug reports, and issues relating to installing or running this software at Purdue University  to the primary maintainer:
* [Michael Hayashi](mailto:mhayashi@purdue.edu?subject=Inquiry%20for%20gds2Para), Graduate Research Assistant, School of ECE, Purdue University

## Installation and Usage
1. Clone this respository into your working location with `git clone git@github.com:purdue-onchip/gds2Para.git`
  * If cloning fails, the usual reason is that a secure connection to the server could not be made. It is important to clone to be able to use `git` to receive updates. These steps will help get your going using `git`:
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
  * Move the zip archive to the working location
  * Unzip the archive with `unzip gds2Para-master.zip`
  * Rename the recently created directory from "gds2Para-master" to "gds2Para"
2. Clone the [Limbo repository](https://github.com/limbo018/Limbo) into your working location with `git clone https://github.com/limbo018/Limbo.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Limbo repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip Limbo-master.zip`
  * Rename the recently created directory from "Limbo-master" to "Limbo"
3. Clone the [Parser-SPEF repository](https://github.com/OpenTimer/Parser-SPEF) into your working location with `git clone https://github.com/OpenTimer/Parser-SPEF.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Parser-SPEF repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip Parser-SPEF-master.zip`
  * Rename the recently created directory from "Parser-SPEF-master" to "Parser-SPEF"
4. Clone the [eigen-git-mirror repository](https://github.com/eigenteam/eigen-git-mirror) into your working location with `git clone https://github.com/eigenteam/eigen-git-mirror.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the eigen-git-mirror repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip eigen-git-mirror-master.zip`
  * Rename the recently created directory from "eigen-git-mirror-master" to "eigen-git-mirror"
5. Ensure that you are the owner of the files that were downloaded with full read/write/execute permissions
  * Change ownership with `chown -R <username>:<username> <directory>`, where \<username> is your Purdue username and \<directory> is each of "gds2Para/", "Limbo/", "Parser-SPEF/", or "eigen-git-mirror/"
  * Change permissions with `chmod -R 744 <directory>`, where \<directory> is each of "gds2Para/", "Limbo/", "Parser-SPEF/", or "eigen-git-mirror/"
6. Ensure that GNU bison is installed on the machine by seeing if there is an output to the terminal with `which bison`
7. Modify run commands files depending on the shell indicated by `echo $SHELL`, substituting \<absolute path to Limbo directory> for a valid path when it appears
  * For "bash", edit the file ".bashrc" in your home directory by appending the following:
```bash
# Skip rest of file if not interactive
if [ -z "$PS1" ]; then
    return
fi

# Improve Compilers
alias cc="gcc"
export CC=/opt/gcc/7.1.0/bin/gcc
export CXX=/opt/gcc/7.1.0/bin/gcc
export FC=/opt/gcc/7.1.0/bin/gfortran
export PATH=/opt/gcc/7.1.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/gcc/7.1.0/lib64:/opt/gcc/7.1.0/lib:$LD_LIBRARY_PATH

# Improve OpenMPI Version
export PATH=/opt/mpich2-gcc/1.4.1p1/bin:${PATH}
export LD_LIBRARY_PATH=/opt/mpich2-gcc/1.4.1p1/lib64:${LD_LIBRARY_PATH}
export OMPI_CC=/opt/gcc/7.1.0/bin/gcc
export OMPI_CXX=/opt/gcc/7.1.0/bin/gcc
export OMPI_FC=/opt/gcc/7.1.0/bin/gfortran

# Improve git Version
export PATH=/opt/git/2.18.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/git/2.18.0/lib64:/opt/git/2.18.0/libexec:$LD_LIBRARY_PATH

# Environment Variables for GDSII File Handling
export BOOST_DIR="/opt/boost/1.57.0"
export FLEX_DIR="/usr/bin/flex"
export LIMBO_DIR="<absolute path to Limbo directory>"
export MKL_DIR="/opt/intel/current/mkl"
```
  * For "tcsh", edit the file ".cshrc" in your home directory by appending the following:
```tcsh
# Skip Rest of File if Not Interactive
if (! $?prompt) then
    exit 0
endif

# Improve gcc Version
alias cc 'gcc'
setenv CC /opt/gcc/7.1.0/bin/gcc
setenv CXX /opt/gcc/7.1.0/bin/gcc
setenv FC /opt/gcc/7.1.0/bin/gfortran
setenv PATH /opt/gcc/7.1.0/bin:${PATH}
setenv LD_LIBRARY_PATH /opt/gcc/7.1.0/lib64:/opt/gcc/7.1.0/lib:${LD_LIBRARY_PATH}

# Improve OpenMPI Version
setenv PATH /opt/mpich2-gcc/1.4.1p1/bin:${PATH}
setenv LD_LIBRARY_PATH /opt/mpich2-gcc/1.4.1p1/lib64:${LD_LIBRARY_PATH}
setenv OMPI_CC /opt/gcc/7.1.0/bin/gcc
setenv OMPI_CXX /opt/gcc/7.1.0/bin/gcc
setenv OMPI_FC /opt/gcc/7.1.0/bin/gfortran

# Improve git Version
setenv PATH /opt/git/2.18.0/bin:${PATH}
setenv LD_LIBRARY_PATH /opt/git/2.18.0/lib64:/opt/git/2.18.0/libexec:${LD_LIBRARY_PATH}

# Environment Variables for GDSII File Handling
setenv BOOST_DIR /opt/boost/1.57.0
setenv FLEX_DIR /usr/bin/flex
setenv LIMBO_DIR <absolute path to Limbo directory>
setenv MKL_DIR /opt/intel/current/mkl
```
8. Exit the shell and terminate the connection before logging back in
9. Ensure that the run command files were properly loaded by running `echo $LIMBO_DIR`
  * If nothing shows up for "bash" users, run `cp .bashrc .bash_profile` in the home directory, exit the shell, log back in, and try again
  * If nothing shows up for "tcsh" users, run `cp .cshrc .tcshrc` in the home directory, exit the shell, log back in, and try again
  * For all other errors, contact the primary maintainer
10. Enter the directory "Limbo/limbo/parsers/gdsii/stream" from the working location
11. Run `make` in this directory to build Limbo libraries
12. Return to working directory and ensure that a new library file named **libgdsparser.a** exists by running `ls -lh Limbo/lib`
13. Follow the instructions in [HYPRE Setup](#HYPRE-Setup) to prepare multigrid methods from HYPRE
  * It is necessary to make all changes to ".bashrc" and ".cshrc" in the home directory before following these instructions
  * Continue following these steps once HYPRE is installed
14. Enter the directory "gds2Para" from the working location
15. (Optional) Add GDSII file besides "nand2.gds", "SDFFRS_X2.gds", and "4004.gds" to working location
  * Some usage modes of this software may also require a simulation input file in a similar format to "nand2.sim_input", "SDFFRS_X2.sim_input", and "4004.sim_input" (the file name before the suffix must match that of the GDSII file before ".gds")
  * Files created from the software will be added directly to the present working directly, which is "gds2Para" within the working location for these instructions
16. Run `make` in shell to compile executable **LayoutAnalyzer**
17. Run `LayoutAnalyzer --help` to get a list of available modes the executable supports
18. Run `LayoutAnalyser -r nand2.gds` in shell to produce terminal output describing the GDSII file ("nand2.gds" may be replaced with any GDSII file available)
19. Perform a complete parameter extraction by running `LayoutAnalyzer -s SDFFRS_X2.gds SDFFRS_X2.sim_input SDFFRS_X2.cir` to read in the design and simulation input file, do all analysis, and return the results in a SPICE subcircuit

## HYPRE Setup
1. Clone the [HYPRE repository](https://github.com/LLNL/hypre) into your working location with `git clone git@github.com:LLNL/hypre.git`
  * This requires the setup of an SSH key for a secure connection
  * Setup will follow the instructions at "hyper/INSTALL" using information known to affect Purdue ECN users
2. Enter the directory "hypre/src/cmbuild" from the working location
  * Configuration options may be viewed using `../configure --help`
<!---3. Run the command `./configure --enable-shared --enable-complex --with-MPI=/opt/mpich2-gcc/1.4.1p1 --with-MPI-include=/opt/mpich2-gcc/1.4.1p1/include/ --with-MPI-libs="" --with-MPI-lib-dirs=/opt/mpich2-gcc/1.4.1p1/lib64/ CC=/opt/gcc/7.1.0/bin/gcc CFLAGS="-g" CXX=/opt/gcc/7.1.0/bin/gcc CXXFLAGS="-g" FC=/opt/gcc/7.1.0/bin/gfortran`--->
3. Run `cmake --version` and ensure that the version of CMake installed is at least 2.8.8
4. Run `cmake -LA ..` and ensure that CMake found the correct versions of the compilers
  * CMake is rather temperamental and tends to overlook common methods of specifying different compilers system-wide
  * Check ".bashrc" or ".cshrc" in your home directory to see in the environment variables for the compilers are properly set
  * For all other issues, contact the primary maintainer
5. Run `cmake ..` to produce a makefile
6. Run `make install` from within the same location ("hypre/src/cmbuild" in the working location) to install header files
7. Return to working directory and ensure that a new library file named **HYPRE.h** exists by running `ls -lh hypre/src/hypre/include`
