#!/bin/csh

# Get's the latest version of Fastjet that Maurizio has in his area,
# and do some preparation.
echo "Getting Fastjet from Maurizio area."
#scp /afs/cern.ch/user/m/mpierini/public/fastjet-2.4.1.tar.gz .
tar -xzf fastjet-2.4.1.tar.gz 
mkdir FASTJET

setenv CXXFLAGS "-m64"
setenv FFLAGS "-m64"
setenv CPPFLAGS "-m64"
setenv CFLAGS "-m64"

set DIR = $PWD/FASTJET
# Compile and install Fastjet.
cd fastjet-2.4.1
./configure --prefix=$DIR
make 
make check
make install

# Come back to the original directory, and clean up.
cd ../
\rm -r fastjet-2.4.1
\rm fastjet-2.4.1.tar.gz

set version = `$PWD/FASTJET/bin/fastjet-config --version`
echo "*******************************************************************"
echo "Fastjet version installed in :$PWD/FASTJET : $version"
echo "*******************************************************************"

#source setup_root_64.csh
