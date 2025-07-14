## VIV Numerical Codes


This repo provide two parts VIV code:

- Flexible pipe VIV of linearly sheared flow case
  - fluids folder contains strips openfoam setting with of8
  - when time step is getting smaller, the precision should be changed to a higher number
  - running `main.m` for simulation
  - `main_continue.m` for resuming calculation from breakpoint
  - reference paper:
      - Fu X, Fu S, Niu Z, et al. A validated fluid-structure interaction simulation model for vortex-induced vibration of a flexible pipe in steady flow[J]. Marine Structures, 2025.



- Bluff cylinder 2DOF case
    - this is the single 2DOF setting for the OpenFOAM-8 to achive 2DOF VIV simulation based on `dynamicMotionSolverFvMesh`
    - this repo is the freely vibrating cylinder case
    - foced motion can be easily achived based on the setting.
    - reference papers:
      - Fu X, Fu S, Han Z, et al. Numerical simulations of 2-DOF vortex-induced vibration of a circular cylinder in two and three dimensions: A comparison study[J]. Journal of Ocean Engineering and Science, 2023.
      - Fu X, Fu S, Zhang M, et al. Frequency capture phenomenon in tandem cylinders with different diameters undergoing flow-induced vibration[J]. Physics of Fluids, 2022, 34(8).
  


