## Installation for PrimerTK and third party software

The easiest way to install PrimerTK is to clone the github repository. This is because the repository contains template scripts for input files, a helper script to run in silico pcr, and a helper script to install primer3 and in silico pcr.

To clone the repo in your current working directory:

```
# https
git clone https://github.com/stjude/PrimerTK.git

# git
git clone git@github.com:stjude/PrimerTK.git
```

Once it is cloned, you can change to the directory and see what is inside:

```
cd PrimerTK
ls
```

This will show the directories `cwl`, `scripts`, `src/primer_tk`, `input_templates`, and `test`. If you read the [introduction](introduction.md), you will know that cwl is a workflow language that will operate the entire pipeline for you after all programs have been installed and are in your `PATH` environment variable.

`input_templates` is just a helper directory for you to see what input files should look like for the program.

To install the python source code to your local python3 using the git copy (assuming you are still in the PrimerTK directory):

```
# to your user python3
pip3 install . --user

# to root python3
pip3 install .
```

This will install the PrimerTK python code. To check for a proper installation, you can type `primer_tk -h` in the command line which will return:

```
usage: primer_tk [-h] [-v]
                 {iterator,iterator_sv,pre,pre_sv,post,post_sv,tabix} ...

positional arguments:
  {iterator,iterator_sv,pre,pre_sv,post,post_sv,tabix}
                        Actions
    iterator            Iterator subparser
    iterator_sv         iterator_sv subparser
    pre                 Preprocessing for snv/indel
    pre_sv              Preprocessing for SV's
    post                Parses info from pcr_output
    post_sv             Parses info from pcr_output
    tabix               Tabix subparser

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

Upon proper installation.

To install third party software, change into the scripts directory and open `installer.sh` with your favorite text editor and remove all commented lines in between text boxes so the file looks like this:

```bash
#
# Third party tools Primer3 (TODO: Link), available under GNU GPLv2, and In Silico PCR
# (TODO: Link) are needed for PrimerTK to work properly. We do not include those third
# party tools alongside the PrimerTK software. However, the following commands were
# used to install the software as of the time of PrimerTK's release.

# Uncomment all code in between boxes to install

###########
# Primer3 #
###########

#############################################################
# Select the installation directory and make sure it exists.#
#############################################################

PRIMER3_INSTALL_DIR=~/bin/
mkdir -p $PRIMER3_INSTALL_DIR

#####################################################################
# Make a tmp location to download / install the source              #
# we need the whole source to run the program. There are configs in #
# primer3_config that are referenced for thermodynamic params       #
#####################################################################

TMPDIR=$(mktemp -d)
echo $TMPDIR
cd $TMPDIR
git clone https://github.com/primer3-org/primer3.git primer3
cd ./primer3/src
make

####################################
# Copy the appropriate executables #
####################################

cp -rf $TMPDIR/primer3 $PRIMER3_INSTALL_DIR

rm -rf $TMPDIR

#########
# isPCR #
#########

ISPCR_EXECUTABLE_DIR=~/bin/isPcr33
mkdir -p $ISPCR_EXECUTABLE_DIR
MACHTYPE=$(uname -m)
export MACHTYPE
ISPCR33_INSTALL_DIR=~/bin/$MACHTYPE
mkdir -p $ISPCR33_INSTALL_DIR

TMPDIR2=$(mktemp -d)
echo $TMPDIR2
cd $TMPDIR2
wget https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip && unzip isPcr33.zip 1> .unzip.out 2> .unzip.err
mkdir ./isPcrSrc/lib/$MACHTYPE
cd $TMPDIR2/isPcrSrc/lib && make HG_WARN=""
cd $TMPDIR2/isPcrSrc && make HG_WARN=""
cp $ISPCR33_INSTALL_DIR/isPcr $ISPCR_EXECUTABLE_DIR
rm -rf $ISPCR33_INSTALL_DIR
rm -rf $TMPDIR2

#################
# END INSTALLER #
#################
```

primer3 and isPcr will now be in `~/bin/`.

I would recommend adding them to your `PATH` variable as they should not conflict with other programs.

You can do this manually every time you run the program, or you can add it to your `.bashrc`

To do it manually:

```
export PATH=~/bin/primer3/src/:$PATH
export PATH=~/bin/isPcr33/:$PATH
```

Or correspondingly you can add the same two lines of code to your `~/.bashrc` and these will both be added to your `PATH` when you login to your account.

This will complete the installation. Next, visit [inputs](inputs.md) to see the structure of the various input files.
