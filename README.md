# Heinz

## Compilation Instructions

### Dependencies

#### Libraries

* LEMON 1.3.1
* IBM ILOG CPLEX >= 12.0, works at least up to 22.1.0
* OGDF 2015-06-29

### Compiling

#### macOS Specifics

You need to install `cmake`, e. g. using `brew`, and the Command Line Developer Tools, which you already have if you installed `cmake` using `brew`. MacOS works at least up to Monterey 12.4.

#### Linux Specifics

On a Debian-derived Linux you need to install `git`, `cmake`, `make`, `g++`, `unzip`, `curl`, and `ca-certificates`. Debian 9 to 11 are known to work as well as Ubuntu 22.04.

#### General Instructions

This guide shows you how to compile heinz, and its dependencies. Because OGDF has to be an old version, the guide places it into a non-standard location, so that it doesn't interfere with the rest of your system.

Get CPLEX from the IBM web-site. For reference, the _Part Number_ for Linux is M04HKML. For macOS it's M04HMML. If you install it into its default location, it should automatically be detected during the `cmake` invocation of the heinz step.

First, LEMON 1.3.1 needs to be installed:

    curl https://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz | tar --extract --gunzip
    cd lemon-1.3.1
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/heinz-dependencies
    make install
    cd ..

Once LEMON is installed, you can remove the source files, i. e. `rm -rf lemon-1.3.1`.
Next, OGDF needs to be built:

    curl --remote-name https://ogdf.uos.de/wp-content/uploads/2019/04/ogdf-snapshot-2015-06-29.zip
    unzip ogdf-snapshot-2015-06-29.zip
    cd OGDF-snapshot
    cmake .
    make
    cp -r include/* $HOME/heinz-dependencies/include
    cp *.a $HOME/heinz-dependencies/lib/
    cd ..

You can remove the OGDF sources now, i.e., `rm -rf OGDF-snapshot ogdf-snapshot-2015-06-29.zip`.

Get heinz, and compile it:

    git clone https://github.com/ls-cwi/heinz
    cd heinz
    mkdir build
    cd build
    cmake -DLIBLEMON_ROOT=$HOME/heinz-dependencies -DLIBOGDF_ROOT=$HOME/heinz-dependencies ..
    make

In case the auto-detection of LEMON, OGDF, or CPLEX fails, you can manually specify they paths:

    cmake \
    -DLIBLEMON_ROOT=~/lemon \
    -DLIBOGDF_ROOT=~/ogdf \
    -DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
    -DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
    -DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
    -DCONCERT_INC_DIR=~/ILOG/concert/include/ ..

When you are done compiling heinz, you can remove the dependencies, i. e. `rm -rf $HOME/heinz-dependencies`.

## Running heinz

To run heinz on the test instances:

    ./heinz -n ../data/test/NodesPCST.txt -e ../data/test/EdgesPCST.txt

To run heinz on the DIMACS MWCS instances, do:

    ./heinz -stp ../data/DIMACS/mwcs/ACTMOD/HCMV.stp

For the PCST DIMACS instances use:

    ./heinz -stp-pcst ../data/DIMACS/pcst/PCSPG-JMP/K100.2.stp
