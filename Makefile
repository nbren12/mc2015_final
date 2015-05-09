
build/param_setter: build/CMakeCache.txt
	cd build ; make

build/CMakeCache.txt:
	mkdir build && cd build && cmake ../src

