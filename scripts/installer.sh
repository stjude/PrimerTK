#
# Third party tools Primer3 (TODO: Link), available under GNU GPLv2, and In Silico PCR
# (TODO: Link) are needed for PrimerTK to work properly. We do not include those third
# party tools alongside the PrimerTK software. However, the following commands were 
# used to install the software as of the time of PrimerTK's release.


###########
# Primer3 #
###########

# Select the installation directory and make sure it exists.
PRIMER3_INSTALL_DIR=~/bin/primer3
mkdir -p $PRIMER3_INSTALL_DIR

# Make a tmp location to download / install the source
TMPDIR=$(mktemp -d)
echo $TMPDIR
cd $TMPDIR
git clone https://github.com/primer3-org/primer3.git primer3
cd ./primer3/src
make
#make test

# Copy the appropriate executables
cp $TMPDIR/primer3/*.pl $PRIMER3_INSTALL_DIR
cp $TMPDIR/src/primer3* $PRIMER3_INSTALL_DIR
# TODO Any more needed?

# Remove the tmp space
#rm -rf $TMPDIR

#########
# isPCR #
#########

# Select the installation directory and make sure it exists
ISPCR_INSTALL_DIR=~/bin/isPCR
mkdir -p $ISPCR_INSTALL_DIR

# Make a tmp location to download / install the source
MACHTYPE=$(uname -m)
TMPDIR=$(mktemp -d)
echo $TMPDIR
cd $TMPDIR
wget https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip
unzip isPcr33.zip 1> .unzip.out 2> .unzip.err
cd $TMPDIR/isPcrSrc/lib && make HG_WARN=""
cd $TMPDIR/isPcrSrc && make HG_WARN=""

# TODO - why doesn't this compile?


