#!/bin/sh
#
# script runs solver 
#

echo "Starting at: `date`"
currentdir=`pwd`

echo "Setup infrastructure..."
mkdir -p input
mkdir -p output
mkdir -p bin
rm -rf output/*
mkdir -p output/velocity
mkdir -p output/stress

echo "Compiling program..."
cd ../../
#make cleanall
make all
if [ $? -ne 0 ]
then
    echo "make failed. Abort..."
    cd $currentdir
    exit 1
fi
cd $currentdir

echo "Copying binaries..."
cd bin
cp ../../../bin/solver .
cd ..

echo "Storing input parameters..."
cp input/* output/

echo "Execute program..."
./bin/solver

echo "done at: `date`"
