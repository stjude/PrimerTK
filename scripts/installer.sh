#
# Third party tools Primer3 (TODO: Link), available under GNU GPLv2, and In Silico PCR
# (TODO: Link) are needed for PrimerTK to work properly. We do not include those third
# party tools alongside the PrimerTK software. However, the following commands were 
# used to install the software as of the time of PrimerTK's release.


###########
# Primer3 #
###########

# Select the installation directory and make sure it exists.
PRIMER3_INSTALL_DIR=~/bin/
mkdir -p $PRIMER3_INSTALL_DIR

# Make a tmp location to download / install the source
# we need the whole source to run the program. There are configs in
# primer3_config that are referenced for thermodynamic params

TMPDIR=$(mktemp -d)
echo $TMPDIR
cd $TMPDIR
git clone https://github.com/primer3-org/primer3.git primer3
cd ./primer3/src
make
#make test

# Copy the appropriate executables
cp -rf $TMPDIR/primer3 $PRIMER3_INSTALL_DIR

# Remove the tmp space
rm -rf $TMPDIR

#########
# isPCR #
#########

# Select the installation directory and make sure it exists
ISPCR_INSTALL_DIR=~/bin/isPcr33
mkdir -p $ISPCR_INSTALL_DIR
# Make a tmp location to download / install the source
MACHTYPE=$(uname -m)
export MACHTYPE
ISPCR_33=~/bin/$MACHTYPE
mkdir -p $ISPCR_33
export PATH=$PATH:$ISPCR_33

TMPDIR2=$(mktemp -d)
echo $TMPDIR2
cd $TMPDIR2
wget https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip && unzip isPcr33.zip 1> .unzip.out 2> .unzip.err
mkdir ./isPcrSrc/lib/$MACHTYPE
cd $TMPDIR2/isPcrSrc/lib && make HG_WARN=""
cd $TMPDIR2/isPcrSrc && make HG_WARN=""
cp $ISPCR_33/isPcr $ISPCR_INSTALL_DIR
rm -rf $ISPCR_33
rm -rf $TMPDIR2
# TODO - why doesn't this compile?
# there is code in the installer which specifically points to /usr/bin/$MACHTYPE as the place to put executables

