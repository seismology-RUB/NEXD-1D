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
make analytical
cd $currentdir

echo "Copying binaries..."
cd bin
cp ../../../bin/analytical .
cd ..

echo "Storing input parameters..."
cp input/* output/

echo "Execute program..."
./bin/analytical

echo "done at: `date`"
