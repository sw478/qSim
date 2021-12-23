# Usage

### Creating a simulation

1. Instantiate a Sim() object, with the desired number of qubits.
2. Add gates sequentially with add_gate(gate, i_qbits).
    a. gate: Gate object describing the gate used
    b. i_qbits: Ordered set of indices of qubits to be used for the gate
3. Call print_statevector(), print_prob_dist(), print_sim() to print Sim info.

### Creating a custom gate from a Sim

1. Create a Sim with desired gates
2. Call turn_into_gate() to return a Gate object
    a. Optional argument for a gate label
3. Can now use Gate in another Sim

## Project Info

* LSB bit ordering similar to [Qiskit](https://qiskit.org/documentation/tutorials/circuits/3_summary_of_quantum_operations.html#Basis-vector-ordering-in-Qiskit)
* All qubits are initialized to [0, 1].

### Implemented Gates

* [Wikipedia Quantum Logic Gates](https://en.wikipedia.org/wiki/Quantum_logic_gate)
    * Most gates listed here

## Notes

Condensed info from different sources that were used for this project

### Youtube

* [Microsoft's "Quantum Computing for Computer Scientists" Lecture](https://www.youtube.com/watch?v=F_Riqjdh2oM)
    * Good explainer for the "computing" aspect
* [Domain of Science's "Map of Quantum Computing"](https://www.youtube.com/watch?v=-UlxHPIEVqA)
    * Broad overview of the current state of quantum computing, inspiration for direction project can potentially take
* [3b1b's Raising a Matrix to an Exponent](https://www.youtube.com/watch?v=O85OWBJ2ayo&t=1231s)
    * Some related concepts to quantum gate matrices

### Online Quantum Simulators

* [IBM's Quantum Composer](https://quantum-computing.ibm.com/composer)