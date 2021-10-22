# huma-ergodicObj
An implemenetation of Sequential Action Control and Hybrid Shared Control with examples
## Requirements
 - cmake version 2.8.3 or higher
 
       $ sudo apt-get install cmake
 
 - Armadillo, A C++ linear algebra library
 
       $ sudo apt-get install liblapack-dev libblas-dev libboost-dev
       $ sudo apt-get install libarmadillo-dev
 
 - OpenCV
 
       $ sudo apt-get install libopencv-dev python3-opencv

## To build the project

- Change `/home/kzf5356/human-ergodicObj` Line 14 of CMakeLists.txt to the location of the repository on your machine.

- Run cmake to generate the build directory and Makefile

      $ cmake .

- Generate executables in the /build directory

      $ make

*Note that all examples using images as the reference distribution are not set up to build with CMakelists. Images need to be added to this folder and their filenames added to individual scripts.*

## Run one of the examples

    $ ./build/examplename

- Data from the simulation is saved in a text file
