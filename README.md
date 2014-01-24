Compilation instructions
========================

Dependencies
------------

* LEMON 1.3
* ILOG CPLEX (>= 12.0)

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
    
Note: On MaxOS 10.9, comment out the following two lines in `CMakeLists.txt` before `make install`

    #ADD_SUBDIRECTORY(demo) 
    #ADD_SUBDIRECTORY(tools)

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`. Next, Heinz can be compiled:

    mkdir build
    cd build
    cmake ..

In case auto-detection of LEMON or CPLEX fails, do

    cmake \
    -DLIBLEMON_ROOT=~/lemon \
    -DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
    -DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
    -DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
    -DCONCERT_INC_DIR=~/ILOG/concert/include/ ..

Running heinz
=============

To run heinz on the test instances:

    ./heinz -n ../data/test/NodesPCST.txt -e ../data/test/EdgesPCST.txt

Preprocessing can be enabled by specifying `-p`
