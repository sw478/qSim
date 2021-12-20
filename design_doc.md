# Design Doc

## Problem Description

The goal of this project is to create a classical quantum circuit simulator.
This will be mainly used as practice for coding and as a learning experience for quantum computing.

## Usage

In main(), there are calls to different example simulations "sim_*()".
All qubits are initialized to [0, 1].

To create a circuit:
1. Instantiate a Sim() object, with the number of qubits.
2. Add gates in the appropriate order with add_gate().
4. Call print_statevector() or print_prob_dist() to view state.

## Requirements

* Can create basic quantum circuits of at least 5 qubits
* Gates:
    * Single bit:
        * Identity
        * Pauli (X, Y, Z)
        * Hadamard (H)
        * Phase shift: (S, T, S dagger, T dagger)
    * Many bit:
        * Swap and variations:
            * Square root of Swap
            * Imaginary Swap
            * Square root of Imaginary Swap
    * Gates with arguments:
        * Arguments should be evaluated at runtime
        * Phase shift gate
        * Rotational gates: (Rx, Ry, Rz)
        * Ising Coupling gates: (Rxx, Ryy, Rzz)
    * Control gates:
        * Controlled single bit gates
        * Toffoli
        * Controlled Swap
    * Custom gates:
        * Gates created from the current gates in a Sim
        * For black boxes / oracles
        * Can support arguments
        * When printing as a gate within a sim, does not need to show details
            * A printout similar to the Sim printout will be printed separately
* File organization
    * "main.py": Where sims are created
    * "sim.py": Sim class
    * "gates.py": Gate class, includes definitions
    * "test.py": Test cases
* Support custom qubit ordering for many-bit gates
    * Reason: Most of the simulation implementations I've seen online only allow for fixed-shape gates, with the qubits needing to be in order.
    * Qubit indices and their associated rows/cols to swap should be preprocessed at the start of each Sim
* Doesn't need to support simulating measurements
    * Reason: Doesn't add much to the project
* Printing operations:
    * print_sim(): Prints the circuit diagram
        * One qubit per row
        * Prints each gate in a separate column
        * If it's a controlled gate the control bits will be
            differentiated from the controlled bit(s)
    * print_statevector(): Prints the current statevector as an array of complex numbers
    * print_prob_dist(): Prints the current probability distribution as an array of real numbers
    * When calling add_gate(), there is an option to print the full matrix used
* Miscellaneous
    * Sims can be named

## Testing

Currently using several methods to debug:
* Can inspect the statevector, probability distribution, and gate matrices between gates. Test cases are made to compare actual and expected values.
* Can create and compare the same simulations in this program and in an already working environment (ex. Qiskit)
    * Can be used alongside the test cases
* For bit order swapping: can print out the index of qubits to be swapped along with the index of rows/columns in the matrix to be swapped

## Future Work

* Might add ability to use classical bits and operators/gates to enable combining multiple quantum circuits
* Support conjugate transpose of a gate matrix as a modifier