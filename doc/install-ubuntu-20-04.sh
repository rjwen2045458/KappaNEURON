# Install KappaNEURON on Ubuntu 20.04
# Requirements:
# - A miniconda or conda installation
# - The pacakges below
# sudo apt install g++9 gcc-9 cython
#
# Instructions:
# Source this file in the KappaNEURON/doc directory:
# cd doc
# . install-ubuntu-20-04.sh

conda deactivate
conda env remove --name KappaNEURON
conda create -n KappaNEURON python=2.7 pip scipy numpy matplotlib
conda activate KappaNEURON

## Install location
NRNPREFIX=$HOME/nrn/7.4/

wget https://neuron.yale.edu/ftp/neuron/versions/v7.4/nrn-7.4.tar.gz
tar zxvf nrn-7.4.tar.gz
cd nrn-7.4
patch -p1 < ../neuron7-4-cython.patch
patch -p1 < ../neuron7-4-abs.patch

CC=gcc-9 CXX=g++-9 CXXFLAGS=-Wno-narrowing ./configure --prefix=$NRNPREFIX/nrn --without-iv  --with-nrnpython=python2.7 --without-memacs
make -j 8
make install

pip install KappaNEURON
export PYTHONPATH=${HOME}/nrn/7.4/nrn/lib/python
python2.7 -i -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalcium
