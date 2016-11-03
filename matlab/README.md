**Compilation**

1. Compile with: ``mex CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp -msse2 -msse4" LDFLAGS="\$LDFLAGS -fopenmp" findex.cc -I../ -largeArrayDims``

2. Run any matlab program with efanna! 

-----

**Dependencies**

* Matlab 2010b or above.
* Other prerequisites are the same as C++ prerequisites, see <https://github.com/fc731097343/efanna/blob/master/README.md>

-----

**Examples**

* Example programs are provided under folder "samples". 

* For instance, use ``matlab -nodesktop -nosplash -r "run('./samples/example_buildgraph')"`` to try. (Don't forget to provide inputs! Default inputs are placed under ``~/data/sift/``. You may change this path in the \*.m programs)

* Every sample does the same job as C++ programs under ``../samples/`` do

* Look inside samples for more API details. 

-----

**FAQ**

Common errors and their solutions:

* If error ``redeclaration of C++ built-in type ‘char16_t’`` occurs in compilation
    * Try uncommenting lines starts with ``define``  and ``undef`` at the beginning of BOTH ``findex.cc`` and ``handle_wrapper.hpp``.
    * Compile again. Problem solved.
    * It may indicates you are using very old version of matlab and not fitting your g++ version well.

* If error related to ``{mablab root}/sys/os/glnxa64/libstdc++.so.6`` occurs in runtime (typically, ``version 'CXXABI_1.x.x' not found (required by ...)``, or `` version `GLIBCXX_3.x.x' not found (required by ...)``
    * Try ``export LD_PRELOAD=$LD_PRELOAD:/usr/lib/x86_64-linux-gnu/libstdc++.so.6``, or ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/libstdc++.so.6``. You may need to rewrite the command with your own path of ``libstdc++.so.6`` according to the place you setup your gcc compiler's lib.
    * OR, directly substitute ``/usr/local/MATLAB/{your version}/sys/os/glnx64/libstdc++.so.6`` with your compiler's ``libstdc++.so.6``. You may need to rewrite the command with the place you setup your matlab.
    * It indicates your mex compiler automatically links the g++ lib under your matlab directory, which is not compatible with your g++ compiler, during the compilation step. 
