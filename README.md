Kalman Filter
===========
Given a stream of noisy input data, the Kalman filter provides statisically optimal estimations of the states of a system. The implementation here demonstrates a real-time tracking system using Kalman filter. The green cross presents an optimal estiation of the mouse input postion (the red square), while adding measurement error in the range of [-50, 50] pixels.

The graphic interface uses OpenGL libraries.


Tracking a point:
![alt tag](https://raw.githubusercontent.com/yanyanggithub/codesamples/gh-pages/images/tracking_point.gif)


Tracking trajectory:
![alt tag](https://raw.githubusercontent.com/yanyanggithub/codesamples/gh-pages/images/tracking_trajectory.gif)

Compile the source using gcc:
> g++ main.cpp -o main  -lglut

Run demo
> ./main

