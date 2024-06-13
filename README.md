# Boeing 747 Front Landing Gear

This document details a Matlab Simulink simulation of the Boeing 747 front landing gear, covering assumptions, model structure, and load calculations. 
It includes a control system with PID feedback, geometric modeling, and analysis of weight, aerodynamic forces, and drag strut spring. 
The report features mathematical formulas and diagrams, supported by a GitHub repository for further resources.



## Running instructions

### Setting up

Open and run `b747_front_land_gear.m` and `sas_loads.m` in MATLAB. 
Make sure You have `sas_aero.csv` in the same directory. 
This will generate data for lookup tables. 
You should see plots showing the data on the screen. 

Then open Simulink and load `b747_front_land_gear_sim.slx` file. 
Click run button.

### Changing parameters

Control input (target angle) and air speed can be modified in the Simulink project by changing corresponding blocks (a few examples provided).
PID coefficients can be changed in the PID block in `b747_front_land_gear_sim.slx` Simulink file to achieve different control response.
Parameters for the geometry can be changed by editing assumed values of the constants in the `b747_front_land_gear.m` and `sas_loads.m`. 
Also one can load different aerodynamic data in `sas_aero.csv` file.
For advanced modifications the actuator constants in the Actuator subsystem can be changed in the `b747_front_land_gear_sim.slx` Simulink file.

### Results analyses

The most convenient way to view the simulation results and explore model behavior is by using Data inspector. 
There also are other possibilities, like oscilloscope (already connected) or some custom solutions using additional Simulink functionalities.
