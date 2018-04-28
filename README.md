"Senior year project - Discrete Time Delta-Sigma Modulator Design" 

sd_test
=> Basic NTF design for a 1-bit, 2nd-order Delta-Sigma modulator.

sd_test_func
=> Auxiliary function for sd_test. Simulates and plots results for a given NTF design for a 1-bit, 2nd-order Delta-Sigma modulator.

test_demo_2_a_theoretical
=> Simulates and plots results for all possible variations in the Leslie-Singh topology design. Takes a flag variable as input.
=> flag = 0 -> Leslie-Singh topology (ideal integrators).
=> flag = 1 -> Leslie-Singh topology (ideal integrators) with input amplitude variation.
=> flag = 2 -> Leslie-Singh topology (ideal integrators) with OSR variation.
=> flag = 3 -> Leslie-Singh topology (non-ideal integrators) with integrator gain and pole variation.
=> flag = 4 -> Non-ideal integrator characterization.
=> flag = 5 -> 1-bit, 2nd-order Delta-Sigma modulator (ideal and non-ideal integrators).
=> flag = 6 -> Leslie-Singh topology (ideal integrators) with coefficient mismatch.

test_demo_2_a_pso_snr2
=> Implements the PSO algorithm for the estimation of parameters in a Delta-Sigma modulator instantiated as a part of the Leslie-Singh topology.

pso_gfunc2
=> Auxiliary function for test_demo_2_a_pso_snr2. Implements the update and cost function calculation equations in the PSO algorithm.

