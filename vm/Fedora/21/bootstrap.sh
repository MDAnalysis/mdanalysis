#!/usr/bin/env bash
# Packages for MDAnalysis on Fedora 21

echo "provisioning Python development stack..."
yum install -y python-pip python-devel python-setuptools gcc

echo "provisioning scientific Python stack..."
yum install -y numpy scipy python-biopython python-matplotlib python-networkx ipython python-nose 

echo "provisioning NETCDF4 for Python (including hdf5 libraries)..."
yum install -y netcdf4-python

# support
#yum install -y clustal-omega

echo "provisioning MDAnalysis..."
pip install MDAnalysis MDAnalysisTests


 

