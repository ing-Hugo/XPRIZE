## Phase I: Initial Development and Proof of Concept

### Overview
Phase I of the project focuses on laying the foundation for quantum-enhanced protein folding to support sustainable insect-based nutrition. This phase is dedicated to establishing a proof of concept, validating the approach with small-scale simulations, and setting the stage for future scalability. The work began with the recognition of the global need for sustainable protein sources amid a projected population of 10 billion by 2050, where insects like crickets and mealworms offer a low-impact alternative to traditional livestock.

### Key Objectives
1. **Algorithm Development**: Implement and test the Variational Quantum Eigensolver (VQE) algorithm to compute ground state energies of protein Hamiltonians.
2. **Initial Simulations**: Conduct simulations on small quantum systems (4 to 50 qubits) to model folding trajectories of insect proteins, such as cricket hemoglobin.
3. **Validation Framework**: Establish a hybrid quantum-classical approach, integrating VQE with classical molecular dynamics (e.g., OpenMM) for structural validation.
4. **Feasibility Assessment**: Evaluate the practicality of quantum advantage over classical methods and identify technical challenges.

### Current Work
* **VQE Implementation**: We are developing VQE to approximate ground state energies. Initial results show energies of -1.4477, -3.8494, -2.8133, -4.1935, and -12.7209 Hartree for 4, 8, 14, 30, and 50-qubit systems, respectively, using tools like Qiskit 1.4.2 on Google Colab.
* **Hamiltonian Construction**: Defining protein Hamiltonians with \( Z \), \( X_i X_{i+1} \), and \( Y_0 Y_{L-1} \) terms to represent folding interactions, starting with a 4-qubit model.
* **Hybrid Validation**: Testing the 1VII peptide (35-residue) folding with OpenMM, achieving an RMSD of 0.0784 Å, to correlate quantum energies with classical structures.
* **Benchmarking**: Attempting classical validation with Density Matrix Renormalization Group (DMRG),(expected -13.0 Hartree for 50 qubits).

### Challenges in Phase I
* **Technical**: Noise in Noisy Intermediate-Scale Quantum (NISQ) devices and DMRG validation errors require HPC resources (1000-10,000 CPU cores, 1 TB RAM).
* **Data Access**: Lack of empirical insect protein folding data necessitates partnerships (e.g., with GAIN or FAO).
* **Resource Constraints**: Limited quantum hardware access slows algorithm refinement.

### Outcomes and Next Steps
By the end of Phase I (targeting late 2025), we aim to:
* Confirm VQE’s effectiveness with a validated 50-qubit simulation.
* Prepare a scalable framework for the 2026 prototype.
* Identify collaboration needs to address data and computational gaps. This phase sets the groundwork for Phase II, where we will scale to 100-200 qubits and conduct pilot studies.
Prepare a scalable framework for the 2026 prototype.
Identify collaboration needs to address data and computational gaps.
This phase sets the groundwork for Phase II, where we will scale to 100-200 qubits and conduct pilot studies.
