# koopman-op
An implemenetation of data-driven control using the Koopman Operator.
## Requirements
 - cmake version 2.8.3 or higher
 
       $ sudo apt-get install cmake
 
 - Armadillo, A C++ linear algebra library
 
       $ sudo apt-get install liblapack-dev libblas-dev libboost-dev
       $ sudo apt-get install libarmadillo-dev
 
 - OpenCV
 
       $ sudo apt-get install libopencv-dev python3-opencv

## To build the project

- Change `/home/kzf5356/iiwa_ros/koopman-op` Line 14 of CMakeLists.txt to the location of the repository on your machine.

- Run cmake to generate the build directory and Makefile

      $ cmake .

- Generate executables in the /build directory

      $ make

## Run a simple example of a cart-pendulum that builds a Koopman operator model of the system

    $ ./build/errsac_cp

## Run a simple example of a cart-pendulum that builds a Koopman operator and uses it to calculate a SAC controller

    $ ./build/errsac_koopcp

- Data from the simulation is saved in a text file

## Run an example of the active learning controller on a simulated quad-copter

	$ ./build/al_quad
	
- Data from the simulation is saved in a text file. The control inputs are plotted and the body centered gravity vector is plotted in figure 1.

### To set up your own active learning controller
- You can use the basis template provided in this repo. 
- See al-quad.cpp lines 51--64 for an example of how to instantiate the alk class from src/ActLearnK.hpp. 
	