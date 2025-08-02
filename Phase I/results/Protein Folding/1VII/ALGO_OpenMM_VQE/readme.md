
# Hybrid Protein Folding Algorithm Description

The hybrid protein folding algorithm developed in this study integrates classical molecular dynamics (MD) simulations via OpenMM with quantum variational quantum eigensolver (VQE) computations, leveraging the QuantumRings service for quantum simulation. This methodology aims to synergize the structural precision of classical MD with the energy optimization capabilities of quantum computing, targeting the folding of the 1VII peptide (sequence: MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF) as a model system. The algorithm operates as follows:

## Initialization and Logging
- **Purpose**: Establishes a logging framework to record computational progress and outcomes, ensuring traceability.
- **Functionality**: Creates a directory to store log files timestamped for each simulation run.

## Structural Data Acquisition
- **Purpose**: Retrieves the initial structural template for the target molecule.
- **Functionality**: Downloads the Protein Data Bank (PDB) file for the target molecule (e.g., 1VII) from the RCSB database if not locally available.

## Sequence Analysis
- **Purpose**: Extracts and simplifies the amino acid sequence for folding.
- **Functionality**: Extracts the amino acid sequence from the PDB file and constructs a hydrophobic-polar (HP) model. Hydrophobic residues (M, L, F, V, W, I, Y, A) are labeled ’H’, while others are labeled ’P’, simplifying the folding problem into a lattice-based representation.

## Hamiltonian Construction
- **Purpose**: Formulates a quantum Hamiltonian based on the HP sequence.
- **Functionality**: Maps the HP sequence onto a 6x6 lattice (36 qubits). The Hamiltonian includes pairwise interaction terms between hydrophobic residues, represented as Pauli Z operators, with coupling strengths set to -1.0 for adjacent or lattice-neighboring ’H’ pairs, encoding the energetic favorability of hydrophobic collapse.

## Quantum Circuit Design
- **Purpose**: Designs a parameterized quantum circuit (ansatz) to approximate the ground state.
- **Functionality**: Creates a 36-qubit ansatz using RY rotations and CNOT gates to entangle adjacent qubits. This ansatz is optimized to reflect the protein’s folded configuration.

## Quantum Optimization
- **Purpose**: Minimizes the energy expectation value using the VQE algorithm.
- **Functionality**: Utilizes the QuantumRings service with the Simultaneous Perturbation Stochastic Approximation (SPSA) optimizer. The circuit parameters are iteratively adjusted over 15 iterations, with 512 shots per evaluation, to converge on an optimal energy state (e.g., -7.895537967480735 dimensionless units).

## Position Mapping
- **Purpose**: Maps optimized quantum parameters to spatial coordinates.
- **Functionality**: Assigns each residue an (x, y, z) position in angstroms on the 6x6 lattice, providing a coarse-grained structural hypothesis.

## Classical Validation
- **Purpose**: Refines the VQE-derived positions using classical MD.
- **Functionality**: Integrates VQE positions into an OpenMM simulation with the AMBER14 force field and Particle Mesh Ewald (PME) for nonbonded interactions. The system undergoes energy minimization and a 10,000-step MD simulation at 300 K, yielding a refined energy (e.g., -2263.98 kJ/mol) and structural output (villin.pdb).

## Benchmarking
- **Purpose**: Compares the VQE structure against a baseline MD simulation.
- **Functionality**: Conducts a baseline MD simulation on the native PDB structure, producing an energy (e.g., -2197.84 kJ/mol). The VQE structure is compared for energy differences and, where possible, structural alignment via root-mean-square deviation (RMSD) using VMD, though VMD availability may limit this step.

## Result Synthesis
- **Purpose**: Consolidates and compares energy results.
![vmd_1VII_hibrid_vqe](https://github.com/user-attachments/assets/b95dd468-3de3-4024-81c0-2605114b250a)




- **Functionality**: Combines VQE energy, OpenMM-validated VQE energy, and baseline MD energy for a comparative analysis to evaluate the hybrid approach’s efficacy in predicting protein folding.
- **Note**: This hybrid framework demonstrated a VQE runtime of 98.24 seconds, with the VQE
