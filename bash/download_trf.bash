#!/bin/bash -e

mkdir -p ../bin
cd ../bin
rm -f trf409.linux64
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
chmod +x trf409.linux64
ls -lh ../bin/trf409.linux64
