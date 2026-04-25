# Analytic_Solutions_Spinning-Body-Hamilton-Jacobi
Code for the computation of the analytic solutions to the nearly-equatorial motion of a test particle with precessing spin in Kerr spacetime. The solutions are valid up to linear order in the spin for eccentric and homoclinic orbits.
Companion repo of "Piovano" arXiv:2510.09597 (https://arxiv.org/abs/2510.09597 ) and "Piovano" arXiv:2603.04682 (https://arxiv.org/abs/2603.04682 ) 

Mathematica files
- ```AnalyticSpinOrbitsHamiltonJacobi.wl```: the package includes the analytic spin-corrections to nearly equatorial periodic orbits, homoclinic orbits, and bound plunges.
  
Mathematica files - "Periodic_homoclinic orbits" folder
- ```SpinOrbitsHamiltonJacobi.wl```: package for the calculation of the numerical spin-corrections to generic periodic orbits. It is used for sanity checks. Also available in https://github.com/gabriel-andres-piovano/Spinning-Body-Hamilton-Jacobi.
- ```Spin_corrections_to_orbits_near_equatorial_motion.nb```: presents the functions available in the package ```AnalyticSpinOrbitsHamiltonJacobi.wl``` for periodic and homoclinic motion.
- ```SpinningTrajectory.wl```: Mathematica package from the GitHub repo https://github.com/vskoupy/KerrSpinningFluxes of Viktor Skoupy. This package contains some of the code used in PRD 108, 044041 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.044041 ) (see arXiv:https://arxiv.org/abs/2303.16798 ). It is used for numerical checks for periodic orbits with the analytic results of arXiv:2510.09597
  
Mathematica files - "Boud_plunge" folder
- ```Spin_corrections_to_orbits_near_equatorial_motion_Bound_Plunge.nb```: presents the functions available in the package ```AnalyticSpinOrbitsHamiltonJacobi.wl``` for bound plunges.

Author:
- Gabriel Andres Piovano
