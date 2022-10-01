# MOC_for_Rocket_Nozzle_Design
Method of Characteristics for Rocket Nozzle Design in C++

This program calculates the coordinates of the wall contour for designing the minimum length supersonic nozzle based on the Method of Characteristics. 

In the "moc_v2.cpp" file, the user inputs such as half-throat radius of the nozzle, exit mach number, specific heat ratio, and number of characteristics lines are specified and the output xy coordinates are stored in "nozzle_xy.txt" file. This file can then be modified to import into a CAD software for creating the contour.

The graph of method of characteristics is plotted with the help of matplotlib package of python library ("matplotlibcpp.h").
