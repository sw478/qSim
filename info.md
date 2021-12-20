# Usage

## Run demo

Running main.py will run a demo

## Creating a simulation

1. Instantiate a Sim() object, with the desired number of qubits.
2. Add gates sequentially with add_gate().
3. Call print_statevector(), print_prob_dist() to view state.

## Miscellaneous Project Info

* LSB bit ordering similar to [Qiskit](https://qiskit.org/documentation/tutorials/circuits/3_summary_of_quantum_operations.html#Basis-vector-ordering-in-Qiskit)
* All qubits are initialized to [0, 1].

## Implemented Gates

* [Wikipedia Quantum Logic Gates](https://en.wikipedia.org/wiki/Quantum_logic_gate)
* [Qiskit's Circuit Library](https://qiskit.org/documentation/apidoc/circuit_library.html)

## Quantum Computing Notes

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