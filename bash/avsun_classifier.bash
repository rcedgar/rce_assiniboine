#!/bin/bash -e

cd "$(dirname "$0")/.."
mkdir -p build-linux big_avsun_classifier
g++ -O3 -std=c++17 -Wall -Wextra -pthread \
	-o build-linux/avsun_classifier \
	avsun_classifier/avsun_classifier.cpp
./build-linux/avsun_classifier \
	--avsun big_avsum_thumbs/avsun_thumbnails \
	--decoy big_avsum_thumbs/decoy_thumbnails \
	--out big_avsun_classifier
