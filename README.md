# Railway Hunting Oscillation Simulator

![Simulation Demo](train_simulation.gif)

**Disclaimer:** This project was created for a university assessment at TU Delft. The code is published here for educational and portfolio purposes. If you are a student, do not copy this code for any assignments. 

## Overview
This project is a numerical simulation of the hunting motion of train bogies, which is a dynamic instability that can lead to derailment. 

Trains rely on conical wheels to self-center on tracks. While this mechanism allows trains to navigate turns, it introduces lateral oscillations (hunting). This Python project simulates a train consisting of $n$ connected bogies using a system of coupled differential equations to model their stability on user-defined track geometries. 

The simulation was utilized to determine optimal track banking angles to minimize flange contact time (the duration wheels grind against the rail), improving safety and reducing wear.

## Key Features
- **Multi-Body Physics:** Simulates $n$ bogies connected by springs and dampers, modeling the complex interaction between carriages.
- **Track Geometry Engine:** Parametric generation of custom track layouts
- **Numerical Integration:** Utilizes the Runge-Kutta 45 (RK45) method via `scipy.integrate` for high-precision solving of equations of motion.
- **Live Visualization:** Real-time animation using `matplotlib` to visualize lateral deviations and train pathing.

## Results
Based on the simulation model, an optimization study was conducted for a quarter-circle turn (Radius = 300m).
- **Findings:** A banking angle of 0.33 mrad was identified as the theoretical optimum.
- **Impact:** This specific angle reduced the flange contact time by 12 seconds compared to a flat track, significantly improving stability.

## Technologies Used
- **Language:** Python 3
- **Computation:** NumPy, SciPy (RK45 Solver)
- **Visualization:** Matplotlib (Animation)

## The Mathematical Model
The simulation is based on Newton's Second Law applied to the lateral ($z$) and yaw ($\phi$) degrees of freedom. The system solves coupled ODEs accounting for:

1.  **Gravitational Restoration:** Forces resulting from track banking angles ($\beta$).
2.  **Spring Coupling:** Interaction forces between adjacent bogies modeled with spring constants ($k$).
3.  **Damping:** Energy dissipation in the suspension system.
4.  **Slipping/Creep:** Modeling the velocity vectors of wheels when lateral oscillation occurs.

## How to Run
1.  Clone the repository:
    ```bash
    git clone [https://github.com/LucasRie/Train-Simulation.git](https://github.com/LucasRie/Train-Simulation.git)
    ```
2.  Install dependencies:
    ```bash
    pip install numpy scipy matplotlib
    ```
3.  Run the simulation:
    ```bash
    python train_simulation.py
    ```
