# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

---

## Architecture
The original code has been refactored to be in the form of a modular library. All the "problem-specific" information, such as the dynamic 
model, the sensor models etc. have been made configurable parts of the UKF class that are passed in to the constructor. 
This can been seen in the first 200 lines of main.cpp. There is no problem-specific information anywhere else in the repository.

The UKF class can utilize any number of sensor models for filtering as long as they are passed in to it during instantiation. 
The dynamic model is also no longer restricted to just CTRV. There is also an option to specify a post processing function for 
both sensor measurements as well as states. This is useful for normalizing angles,putting bounds on values etc.


## *Notes/Caveats*
 
In order to make the filter numerically stable, the maximum turn-rate is bounded to 60 degrees/s and the filter is reset if the
covariance matrix becomes non-positive semi definite at any point.

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.

## Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/c3eb3583-17b2-4d83-abf7-d852ae1b9fff/concepts/4d0420af-0527-4c9f-a5cd-56ee0fe4f09e)
for instructions and the project rubric.
