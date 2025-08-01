# XPRIZE

# Quantum-Enhanced Protein Folding for Sustainable Insect-Based Nutrition

## Overview
This project harnesses quantum computing to transform sustainable nutrition by optimizing the protein folding of insects like crickets and mealworms. As the global population nears 10 billion by 2050, traditional livestock farming burdens resources and emits high greenhouse gases. Insects provide a low-impact alternative, using 80% less land and 100-fold fewer emissions per kilogram of protein. We use quantum algorithms to unlock their nutritional potential for human consumption.

## Project Goals
- Deliver a quantum prototype by 2026 to simulate a 10-50 qubit protein fragment (e.g., cricket hemoglobin).
- Feed 10,000 individuals in a pilot study by 2026.
- Cut livestock emissions by 5-10% by 2050, supporting SDGs 2, 12, and 13.

## Methodology
### Variational Quantum Eigensolver (VQE)
- **Purpose**: Computes ground state energy of protein Hamiltonians.
- **Process**: Employs a parameterized quantum circuit (ansatz) with classical optimization (e.g., SPSA), achieving -1.4477 Hartree for a 4-qubit system.
- **Application**: Models folding for 4, 8, 14, 30, and 50-qubit systems.

### Harrow-Hassidim-Lloyd (HHL) Algorithm
- **Purpose**: Solves linear systems to refine folding dynamics.
- **Process**: Uses quantum phase estimation for \( O(\log n) \) scaling.
- **Application**: Planned for 100-200 qubit systems.

### Hybrid Approach
- Integrates VQE with OpenMM for structural validation, e.g., 1VII peptide with 0.0784 Å RMSD.

## Results
- VQE energies: -1.4477, -3.8494, -2.8133, -4.1935, -12.7209 Hartree for 4, 8, 14, 30, and 50 qubits.


## Getting Started
1. **Prerequisites**: Python 3.x, Qiskit 1.4.2, OpenFermion, OpenMM.
2. **Installation**: Clone repo and run `pip install -r requirements.txt`.
3. **Usage**: Execute `python main.py` for sample VQE simulations.

## Repository Structure
- `/src/`: VQE and HHL algorithm code.
- `/data/`: Hamiltonians and protein data.
- `/results/`: Energy outputs and visualizations.

## Challenges
- Technical: Scaling to 100-200 qubits, NISQ noise.
- Financial: Quantum hardware access.
- Organizational: Team coordination.

## Future Work
- Improve Algorithms with HPC resources.
- Expand to other insect proteins by 2026.

## Contributions
Contributions are welcome! Fork, create a feature branch, and submit a pull request.

## License
[MIT License](LICENSE)

## Contact
visit [https://github.com/ing-Hugo/XPRIZE](https://github.com/ing-Hugo/XPRIZE).
