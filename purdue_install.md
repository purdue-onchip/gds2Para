# Instructions Specific to Purdue ECN Machines
Please follow these instructions for a more complete deployment on Red Hat Enterprise Linux (RHEL) machines at Purdue University maintained by the Engineering Computer Network (ECN). These instructions do **not** apply to any users outside Purdue University.

## Developers
Please direct all questions, bug reports, and issues relating to installing or running this software at Purdue University  to the primary maintainer:
* [Michael Hayashi](mailto:mhayashi@purdue.edu?subject=Inquiry%20for%20gdsii-interface), Graduate Research Assistant, School of ECE, Purdue University

## Installation and Usage
1. Clone this respository into your working location with `git clone git@github.com:purdue-onchip/gdsii-interface.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of this repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip gdsii-interface-master.zip`
2. Clone the [Limbo repository](https://github.com/limbo018/Limbo) into your working location with `git clone https://github.com/limbo018/Limbo.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Limbo repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip Limbo-master.zip`
3. Ensure that you are the owner of the files that were downloaded with full read/write/execute permissions
  * Change ownership with `chown -R <username>:<username> <directory>`, where \<username> is your Purdue username and \<directory> is either "gdsii-interface/" (cloned), "gdsii-interface-master/" (downloaded), "Limbo/" (cloned), or "Limbo-master/" (downloaded)
  * Change permissions with `chmod -R 744 <directory>`, where \<directory> is either "gdsii-interface/" (cloned), "gdsii-interface-master/" (downloaded), "Limbo/" (cloned), or "Limbo-master/" (downloaded)
4. Ensure that GNU bison is installed on the machine by seeing if there is an output to the terminal with `which bison`
5. Modify run commands files depending on the shell indicated by `echo $SHELL`, substituting \<absolute path to Limbo directory> for a valid path when it appears
  * For "bash", edit the file ".bashrc" in your home directory by appending the following:
```bash
# Improve gcc Version
export PATH=/opt/gcc/7.1.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/gcc/7.1.0/lib64:/opt/gcc/7.1.0/lib:$LD_LIBRARY_PATH

# Environment Variables for GDSII File Handling
export BOOST_DIR="/opt/boost/1.57.0"
export FLEX_DIR="/usr/bin/flex"
export LIMBO_DIR="<absolute path to Limbo directory>"
```
  * For "tcsh", edit the file ".cshrc" in your home directory by appending the following:
```tcsh
# Improve gcc Version
setenv PATH /opt/gcc/7.1.0/bin:${PATH}
setenv LD_LIBRARY_PATH /opt/gcc/7.1.0/lib64:/opt/gcc/7.1.0/lib:${LD_LIBRARY_PATH}

# Environment Variables for GDSII File Handling
setenv BOOST_DIR /opt/boost/1.57.0
setenv FLEX_DIR /usr/bin/flex
setenv LIMBO_DIR <absolute path to Limbo directory>
```
6. Exit the shell and terminate the connection before logging back in
7. Enter "Limbo/limbo/parsers/gdsii/stream" (cloned) or "Limbo-master/limbo/parsers/gdsii/stream" (downloaded) directory in working location
8. Run `make CXX=g++ CC=gcc FC=gfortran` in this directory
9. Return to working directory and ensure that a new library file named "libgdsparser.a" exists by running `ls -lh Limbo/lib`
10. Enter "gdsii-interface" (cloned) or "gdsii-interface-master" (downloaded) directory from working location
11. (Optional) Add GDSII file besides "nand2.gds" to working location
12. Run `make explicit` in shell to compile executable **Test\_explicit**
13. Run `./Test_explicit nand2.gds` in shell to produce terminal output describing the GDSII file ("nand2.gds" may be replaced with your own GDSII file if provided earlier)
