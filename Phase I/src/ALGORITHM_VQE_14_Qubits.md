# Pseudo-code for VQE Simulation with Quantum Rings (14 Qubits)

## Overview
This pseudo-code outlines a Variational Quantum Eigensolver (VQE) simulation for a 14-qubit system using a generic quantum backend. It computes an estimated ground-state energy of a Hamiltonian, compares it with a rough lower bound, and visualizes convergence. The code integrates Qiskit and matplotlib, employing an EfficientSU2 ansatz with adaptive optimization.

## Complete Pseudo-code
```plaintext

// Define 14-qubit Hamiltonian
SET hamiltonian_14 = SparsePauliOp FROM LIST [
    ("ZIIIIIIIIIIIII", 0.5), ("IZIIIIIIIIIIII", 0.45), ("IIZIIIIIIIIIII", 0.4),
    ("IIIZIIIIIIIIII", 0.35), ("IIIIZIIIIIIIII", 0.3), ("IIIIIZIIIIIIII", 0.25),
    ("IIIIIIZIIIIIII", 0.2), ("IIIIIIIZIIIIII", 0.15), ("IIIIIIIIZIIIII", 0.1),
    ("IIIIIIIIIZIIII", 0.08), ("IIIIIIIIIIZIII", 0.06), ("IIIIIIIIIIIZII", 0.04),
    ("IIIIIIIIIIIIZI", 0.03), ("IIIIIIIIIIIIIZ", 0.02),
    ("XXIIIIIIIIIIII", 0.2), ("IXXIIIIIIIIIII", 0.18), ("IIXXIIIIIIIIII", 0.16),
    ("IIIXXIIIIIIIII", 0.14), ("IIIIXXIIIIIIII", 0.12), ("YIIIIIIIIIIIIY", 0.1)
]
PRINT "Hamiltonian (14 qubits): " + hamiltonian_14

// Placeholder for Exact Ground State
SET exact_ground_state = -SUM([0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.08, 0.06, 0.04, 0.03, 0.02])
PRINT "Estimated ground-state energy (lower bound): " + exact_ground_state + " Hartree"

SET num_qubits = 14
SET ansatz_14 = EfficientSU2 WITH num_qubits, reps 3, entanglement 'linear'
SET ansatz_transpiled_14 = transpile(ansatz_14, optimization_level 1)

// Custom SPSA Optimizer
DEFINE CLASS CustomSPSA
    DEFINE CONSTRUCTOR(maxiter=3000, learning_rate=0.09, perturbation=0.05)
        SET self.maxiter = maxiter
        SET self.initial_lr = learning_rate
        SET self.initial_pert = perturbation

    DEFINE FUNCTION step(objective_function, params, iteration)
        SET lr = self.initial_lr / (1 + iteration // 100)
        SET pert = self.initial_pert / (1 + iteration // 100)
        SET delta = RANDOM_CHOICE([-1, 1], size params.shape) * pert
        SET loss_plus = objective_function(params + delta)
        SET loss_minus = objective_function(params - delta)
        SET gradient = (loss_plus - loss_minus) / (2 * delta)
        RETURN params - lr * gradient

SET optimizer = CustomSPSA WITH maxiter 3000
SET estimator = Estimator WITH backend
SET estimator.options.default_shots = 4096

// VQE Optimization
DEFINE CUSTOM_VQE(estimator, ansatz, hamiltonian, optimizer)
    SET params = RANDOM(0, 1, ansatz.num_parameters)
    INITIALIZE energies AS EMPTY LIST
    PRINT "Running VQE with " + optimizer.maxiter + " iterations..."
    FOR i FROM 0 TO optimizer.maxiter - 1
        SET pub = (ansatz, hamiltonian, params)
        SET job = estimator.run([pub])
        PRINT "\nIteration " + (i + 1) + "/" + optimizer.maxiter
        SET result = job.result()
        SET energy = result[0].data.evs
        IF energy IS LIST OR np.ndarray
            SET energy = energy[0]
        APPEND energy TO energies
        PRINT "Energy: " + energy + " Hartree"

        DEFINE FUNCTION objective_function(p)
            SET pub = (ansatz, hamiltonian, p)
            SET job = estimator.run([pub])
            SET result = job.result()
            SET energy = result[0].data.evs
            RETURN energy[0] IF energy IS LIST OR np.ndarray ELSE energy
        SET params = optimizer.step(objective_function, params, i)

    SET final_energy = LAST(energies) IF energies IS NOT EMPTY ELSE NONE
    RETURN final_energy, params, energies

PRINT "\nStarting VQE simulation..."
SET final_energy, optimal_params, energy_history = CUSTOM_VQE(estimator, ansatz_transpiled_14, hamiltonian_14, optimizer)

IF final_energy IS NOT NONE
    PRINT "\nGround-state energy: " + final_energy + " Hartree"
    PRINT "Error from exact (estimated): " + (exact_ground_state - final_energy) + " Hartree"
    PRINT "Optimal parameters: " + optimal_params
    PLOT range(1, LENGTH(energy_history) + 1) VS energy_history WITH marker 'o'
    ADD HORIZONTAL_LINE AT exact_ground_state WITH color 'r', linestyle '--', label 'Estimated: ' + exact_ground_state
    SET x-label AS "Iteration"
    SET y-label AS "Energy (Hartree)"
    SET title AS "VQE Convergence (14 Qubits, EfficientSU2 reps=3)"
    ADD legend
    ADD grid
    SHOW plot

PRINT "VQE simulation completed!"
