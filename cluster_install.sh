cd ~

mkdir bin

wget https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip && unzip isPcr33.zip

MACHTYPE=x86_64
export MACHTYPE

mkdir ~/bin/$MACHTYPE

PATH=$PATH:~/bin/$MACHTYPE
export $PATH

mkdir ./isPcrSrc/lib/$MACHTYPE

cd ./isPcrSrc/lib && make HG_WARN=""

wait

cd ../ && make HG_WARN=""

cd ..

cp ~/bin/$MACHTYPE/isPcr .

rm -rf ~/bin/$MACHTYPE

rm -rf ./isPcrSrc/

rm -rf ./isPcr33.zip

