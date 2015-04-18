Compilation instructions
========================

Dependencies
------------

* LEMON 1.3
* ILOG CPLEX (>= 12.0)
* OGDF (v. 2012.07)

Compiling
---------

Get heinz from github:

    git clone <HTTPS clone URL (see on the right side of this page)>


First, LEMON 1.3 needs to be installed:

    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz
    tar xvzf lemon-1.3.tar.gz
    cd lemon-1.3
    cmake -DCMAKE_INSTALL_PREFIX=~/lemon
    make install
    
Note: On Mac OS 10.9, comment out the following two lines and add the code below at line 162 in `CMakeLists.txt` before doing `make install`.

    #ADD_SUBDIRECTORY(demo) 
    #ADD_SUBDIRECTORY(tools)
    
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
    endif()

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`.

Next, OGDF needs to be installed:

    wget http://www.ogdf.net/lib/exe/fetch.php/tech:ogdf.v2012.07.zip
    unzip tech\:ogdf.v2012.07.zip
    cd OGDF
    ./makeMakefile.sh 
    mkdir ~/ogdf
    mkdir ~/ogdf/lib
    mkdir ~/ogdf/include
    make
    cp _release/libOGDF.a ~/ogdf/lib/
    cp -R ogdf ~/ogdf/include/
    
Note: On Mac OS 10.9, before doing `./makeMakefile.sh` edit line 40 in `makeMakefile.config` as follows.

    release = -stdlib=libstdc++ -O2

You can remove the OGDF sources now, i.e., `rm -rf OGDF`. 

Next, Heinz can be compiled:

    mkdir build
    cd build
    cmake ..
    make
    
In case auto-detection of LEMON, OGDF or CPLEX fails, do

    cmake \
    -DLIBLEMON_ROOT=~/lemon \
    -DLIBOGDF_ROOT=~/root \
    -DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
    -DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
    -DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
    -DCONCERT_INC_DIR=~/ILOG/concert/include/ ..

Running heinz
=============

To run heinz on the test instances:

    ./heinz -n ../data/test/NodesPCST.txt -e ../data/test/EdgesPCST.txt

Or do:

    make check
    
To run heinz on the DIMACS MWCS instances, do:

    ./heinz -stp ../data/DIMACS/mwcs/ACTMOD/HCMV.stp
    
For the PCST DIMACS instances use:

    ./heinz -stp-pcst ../data/DIMACS/pcst/PCSPG-JMP/K100.2.stp