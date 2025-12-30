![Simulation Demo](train_simulation.gif)

Disclaimer: This project was created for a university assessment. The code is published here for enducational and portfolio purposes. If you are a student, do not copy this code for any assignments. 

Overview:
This project is a numerical simulation of the hunting motion of train bogies, which assess the stability of railway bogies to determine if the train will derail under certain circumstances. Trains rely on the conical heels to self-center itself on the tracks. This is the main mechanism that allows for trains to take turns, however it comes with the drawback of a dynamic instability known as hunting motion. This python project provides a visualistion of the train consisting of n connected bogies using a systme of coupled differential equations to simulate their movement through a track defined by the user. A further development on this project used this simulation to determine to optimal track banking angles to minimise contact time between the flanges on the wheels and the rails. 

Features:
Multi-Body Physics: Simulates $n$ bogies connected by springs and dampers.
Track Geometry Engine: Parametric generation of track layouts.
Numerical Integration: Uses Runge-Kutta 45 method to solve the equations of motion.
Live Visualization: Real-time animation using matplotlib showing the train's path and lateral deviations.

The mathematical model: The simulation is based on Newton's Second Law applied to the lateral (z) and yaw($\psi$) degrees of freedom. The key forces modelled are: 
1) Gravitational Restoration: a result of track banking.
2) Spring Coupleing: between adjacent bogies.
3) Damping: energy dissipation in the suspension.
4) Creep Forces: interactions between the wheel and rail during slip.
