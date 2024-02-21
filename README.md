This repository contains all FEniCS and Firedrake codes used for our work on reference configuration computations of nonlinear elasticity, found here [todo: add link]

In it you will find: 

## Section 4: Self-intersection of the stress-free state

For this section, we used the ´firedrake/´ folder. 

- For generating the meshes, check the ´memo.txt´ file.
- To generate Figures 4 and 5, use ´inverse-displacement.py´
- To perform the preconditioners comparison, use ´inverse-preconditioning.py´. Here, the third parameter can be either ´hypre´ for AMG or ´mg´ for GDSW.

## Section 6: Numerical tests

For this section, we used the ´fenics/´ folder. 

- To run the Numerical Convergence tests, use the files ´add filenames´.
- To run the robustness tests, use ´run-tests-robustness.sh´. This will spam all tasks through ´tsp´ so that they run in the background. Analogously, use ´run-tests-robustness.sh´ for the other tests.
- To obtain Figure 7, use ´run-tests-display.sh´.
- To compute the performance and CPU time comparisons, use ´run-tests-id-comparison.sh´.
- All tests on the four chambers heart were performed in lifex.





