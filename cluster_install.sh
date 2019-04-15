cd ~

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

module load gcc/6.3.0
module load perl/5.10.1
cd ~

git clone https://github.com/primer3-org/primer3.git primer3

cd ./primer3/src && make
wait
make test

# to add primer3 to your path, do the following:
# cd ~
# emacs .bashrc
# PATH=$PATH:~/primer3/src/
# export $PATH
# upon your next login you will be able to use primer3 from the command line

wait

module load python/3.6.1
pip3.6 install --user 'pandas==0.23.4'
wait
pip3.6 install --user 'numpy==1.16.0'
wait
pip3.6 install --user 'biopython==1.72'
wait
pip3.6 install --user 'cwl-runner'
wait
pip3.6 install --user 'pysam==0.15.2'
wait
pip3.6 install --user 'pyfaidx'

echo '### PrimerPlex paths' >> ~/.bashrc
echo 'module load python/3.6.1' >> ~/.bashrc
echo 'module load gcc/6.3.0' >> ~/.bashrc
echo 'export PATH=$PATH:$(pwd)' >> ~/.bashrc
echo 'export PATH=$PATH:~/.local/bin/' >> ~/.bashrc
echo 'export PATH=$PATH:$(pwd)/primer3/src/' >> ~/.bashrc

