1. Compile with: ``mex CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp -msse2 -msse4" LDFLAGS="\$LDFLAGS -fopenmp" findex.cc -I../ -largeArrayDims``

2. Run any matlab program with efanna! 

-----

**Examples**

- Example programs are provided under folder "samples". For instance, use ``matlab -nodesktop -nosplash -r "run('./samples/example_kdtreeub')"`` to try. (Don't forget to provide inputs! Default inputs are placed under ``~/data/sift/``. You may change this path in the \*.m programs)

- Every sample does the same job as C++ programs under ``../samples/`` do

- Look inside samples for more API details. 

-----

**FAQ**

For those using lower version of matlab and not fitting your gcc version well:

- If error ``redeclaration of C++ built-in type ‘char16_t’`` occurs in compilation, try uncommenting lines starts with ``define``  and ``undef`` at the beginning of``findex.cc`` and ``handle_wrapper.hpp``, then compile again.

- If error ``version 'CXXABI_1.x.x' not found (required by ...)`` occurs in runtime, it indicates the gcc lib version under your matlab directory may be lower than your gcc compiler version in the compilation step. Try exporting your compiler's gcc lib path (which should have ``libstdc++.so.6`` inside) to environment variable ``$LD\_LIBRARY\_PATH``, OR, substitute the gcc lib file in matlab folder, which most probably is ``/usr/local/MATLAB/{your version}/bin/sys/os/glnx64/libstdc++.so.6``, with your compiler's ``libstdc++.so.6``.
