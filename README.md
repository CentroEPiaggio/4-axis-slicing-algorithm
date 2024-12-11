# 4-axis-slicing-algorithm

# Overview
Matlab codes for cylindrical slicing and printing planning on rotating spindle.
These codes allow you to slice an object by dividing it into cylindrical layers defined around a mandrel, rather than into planar layers as in traditional slicers.

# Usage
Like all common slicing software, this one also requires an .STL file of the geometry you want to slice as input.
The code provided considers that the mandrel is oriented with its axis coinciding with the X axis. For this reason, consider rotating your .STL file by orienting it consistently with this axis before loading it into Matlab and slicing.
The zip folder contains all the functions for the correct functioning of the Main Code (Slicer_main). Run the code one section at a time following the instructions and choosing the options contained in the script.
