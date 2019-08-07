#!/usr/bin/env bash

echo "
####################################################
### Primer3 and In Silico PCR License Agreements ###
####################################################

Primer3 is licensed under GNU GPL v2.0. To see in detail with this license entails please visit:
https://github.com/primer3-org/primer3/blob/master/LICENSE. 
This message will serve as a user agreement to download primer3 as a completely independent software to PrimerTK. 
PrimerTK does not serve to distribute primer3 or any of its components. By typing 'yes', you, the user, agree to the
terms and conditions of the GPL2 license and are agreeing to download primer3 independently of PrimerTK. You can
also refuse to download primer3 by entering no.

In Silico PCR [isPCR] is freely available for non-profit/educational use only without a license. 
By typing 'yes' you are agreeing that you are a non-profit/educational institution and do not intend to use this
commercially. We release rights with this installer to the user installing the software and take no responsibility
of actions after install. 

Edit the code and uncomment to install. By uncommenting you are taking ownership of the install.

"


#mydir=$(pwd)
#git clone https://github.com/primer3-org/primer3.git primer3
#cd ./primer3/src && make
#make test
#cd $mydir
#wget https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip && unzip isPcr33.zip
#MACHTYPE=x86_64
#export MACHTYPE
#mkdir ~/bin
#mkdir ~/bin/$MACHTYPE
#export PATH=$PATH:~/bin/$MACHTYPE
#mkdir ./isPcrSrc/lib/$MACHTYPE
#cd ./isPcrSrc/lib && make HG_WARN=""
#cd ../ && make HG_WARN=""
#cd ../
#cp ~/bin/$MACHTYPE/isPcr .
#rm -rf ~/bin/$MACHTYPE
#rm -rf ./isPcrSrc/
#rm -rf ./isPcr33.zip
