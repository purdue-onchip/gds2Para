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
2. Clone the [Limbo repository](https://github.com/limbo018/Limbo) into your working location with `git clone https://github.com/limbo018/Limbo.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Limbo repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip Limbo-master.zip`
3. Clone the [Parser-SPEF repository](https://github.com/OpenTimer/Parser-SPEF) into your working location with `git clone https://github.com/OpenTimer/Parser-SPEF.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the Parser-SPEF repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip Parser-SPEF-master.zip`
4. Clone the [eigen-git-mirror repository](https://github.com/eigenteam/eigen-git-mirror) into your working location with `git clone https://github.com/eigenteam/eigen-git-mirror.git`
  * If cloning fails, try to download the files by clicking the green "Clone or download" button on the Code tab of the eigen-git-mirror repository followed by "Download ZIP"
  * Move the zip archive to the working location
  * Unzip the archive with `unzip eigen-git-mirror-master.zip`
5. Ensure that you are the owner of the files that were downloaded with full read/write/execute permissions
  * Change ownership with `chown -R <username>:<username> <directory>`, where \<username> is your Purdue username and \<directory> is either "gds2Para/" (cloned) or "gds2Para-master/" (downloaded), "Limbo/" (cloned) or "Limbo-master/" (downloaded), "Parser-SPEF/" (cloned) or "Parser-SPEF-master/" (downloaded), or "eigen-git-mirror/" (cloned) or "eigen-git-mirror-master/" (downloaded)
  * Change permissions with `chmod -R 744 <directory>`, where \<directory> is either "gds2Para/" (cloned) or "gds2Para-master/" (downloaded), "Limbo/" (cloned) or "Limbo-master/" (downloaded), "Parser-SPEF/" (cloned) or "Parser-SPEF-master/" (downloaded), or "eigen-git-mirror/" (cloned) or "eigen-git-mirror-master/" (downloaded)
6. Ensure that GNU bison is installed on the machine by seeing if there is an output to the terminal with `which bison`
7. Modify run commands files depending on the shell indicated by `echo $SHELL`, substituting \<absolute path to Limbo directory> for a valid path when it appears
  * For "bash", edit the file ".bashrc" in your home directory by appending the following:
```bash
# Improve gcc Version
alias cc="gcc"
export PATH=/opt/gcc/7.1.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/gcc/7.1.0/lib64:/opt/gcc/7.1.0/lib:$LD_LIBRARY_PATH

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
# Improve gcc Version
alias cc 'gcc'
setenv PATH /opt/gcc/7.1.0/bin:${PATH}
setenv LD_LIBRARY_PATH /opt/gcc/7.1.0/lib64:/opt/gcc/7.1.0/lib:${LD_LIBRARY_PATH}

# Improve git Version                                                                                 
setenv PATH /opt/git/2.18.0/bin:${PATH}                                                               
setenv LD_LIBRARY_PATH /opt/git/2.18.0/lib64:/opt/git/2.18.0/libexec:${LD_LIBRARY_PATH}

# Environment Variables for GDSII File Handling
setenv BOOST_DIR /opt/boost/1.57.0
setenv FLEX_DIR /usr/bin/flex
setenv LIMBO_DIR <absolute path to Limbo directory>
seten MKL_DIR /opt/intel/current/mkl
```
8. Exit the shell and terminate the connection before logging back in
9. Ensure that the run command files were properly loaded by running `echo $LIMBO_DIR`
  * If nothing shows up for "bash" users, run `cp .bashrc .bash_profile` in the home directory, exit the shell, log back in, and try again
  * If nothing shows up for "tcsh" users, run `cp .cshrc .tcshrc` in the home directory, exit the shell, log back in, and try again
  * For all other errors, contact the primary maintainer
10. Enter "Limbo/limbo/parsers/gdsii/stream" (cloned) or "Limbo-master/limbo/parsers/gdsii/stream" (downloaded) directory in working location
11. Run `make CXX=g++ CC=gcc FC=gfortran` in this directory
12. Return to working directory and ensure that a new library file named **libgdsparser.a** exists by running `ls -lh Limbo/lib`
13. Enter "gds2Para" (cloned) or "gds2Para-master" (downloaded) directory from working location
14. (Optional) Add GDSII file besides "nand2.gds", "SDFFRS_X2.gds", and "4004.gds" to working location
15. Run `make` in shell to compile executable **LayoutAnalyzer**
16. Run `LayoutAnalyzer --help` to get a list of available modes the executable supports
17. Run `LayoutAnalyser -r nand2.gds` in shell to produce terminal output describing the GDSII file ("nand2.gds" may be replaced with any GDSII file available)
18. Perform a complete parameter extraction by running `LayoutAnalyzer -s SDFFRS_X2.gds SDFFRS_X2.sim_input SDFFRS_X2.spef` to read in the design and simulation input file, do all analysis, and return the results in SPEF files
