if [ -d "Build" ]; then
	echo "Directory Build already exists. Please delete the whold directory and rerun the script!"
	echo ""
	set -e
fi
mkdir Build
cmake -E chdir Build/ cmake -G "Unix Makefiles" cmake -DCMAKE_BUILD_TYPE=Release ../
