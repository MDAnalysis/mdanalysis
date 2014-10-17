#!/usr/bin/env bash
# Packages for MDAnalysis on Ubuntu 12.04
apt-get update

echo "provisioning Python development stack..."
apt-get install -y build-essential python-dev python-setuptools python-pip

echo "provisioning scientific Python stack..."
apt-get install -y python-numpy python-scipy python-matplotlib python-biopython python-networkx ipython

echo "provisioning NETCDF4 for Python..."
apt-get install -y libhdf5-serial-dev libnetcdf-dev
pip install netCDF4

# support
#apt-get install -y clustalo

echo "provisioning MDAnalysis..."
pip install MDAnalysis MDAnalysisTests


 

