#!/bin/bash
if [ -d "Build" ]; then
	echo "Directory Build already exists. Please delete the whold directory and rerun the script!"
	echo ""
	set -e
fi
mkdir build
cmake -E chdir build/ cmake -G "Unix Makefiles" cmake -DCMAKE_BUILD_TYPE=Release ../
