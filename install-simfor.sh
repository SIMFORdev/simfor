
#!/bin/bash

mkdir simfor
cd simfor
cmake ..
make -j 8
sudo checkinstall -y
cd ..
sudo rm -rfd simfor
