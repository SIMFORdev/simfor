
#!/bin/bash

echo "--- Installing dependencies from debian repositories"
sudo apt update
sudo apt install -y wget make cmake build-essential g++ checkinstall libopenmpi-dev mesa-utils freeglut3-dev doxygen xterm
if [ $? -eq 0 ]; then
	echo "--- Installing dependencies from debian repositories -- OK"
else
  	echo "--- Installing dependencies from debian repositories -- Error"
  	return -1
fi

echo
echo "--- Downloading boost_1_84_0"
cd /tmp/
wget https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.gz
if [ $? -eq 0 ]; then
	echo "--- Downloading boost_1_84_0 -- OK"
else
  	echo "--- Downloading boost_1_84_0 -- Error"
  	rm boost_1_84_0.tar.gz
  	return -1
fi

echo
echo "--- Installing boost_1_84_0"
echo "-- Unpack boost_1_84_0"
tar -xvf boost_1_84_0.tar.gz &> /dev/null
rm boost_1_84_0.tar.gz
cd boost_1_84_0/
echo "-- Configure boost_1_84_0"
./bootstrap.sh
echo "using mpi ;" >> project-config.jam
echo "-- Installing boost_1_84_0"
sudo ./b2 install
if [ $? -eq 0 ]; then
	echo "--- Installing boost_1_84_0 -- OK"
else
  	echo "--- Installing boost_1_84_0 -- Error"
fi

cd /tmp/
sudo rm -frd boost_1_84_0
