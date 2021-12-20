# Design Doc

## Problem Description

The goal of this project is to create a classical quantum circuit simulator.
This will be mainly used as practice for coding and as a learning experience for quantum computing.

## Requirements

* Create quantum circuits similar to Qiskit
* Gates:
    * Single bit:
        * Identity
        * Pauli (X, Y, Z)
        * Hadamard
        * Phase shift: (S, T, S dagger, T dagger)
    * Many bit:
        * Swap and variations (sqrt, imaginary)
    * Gates with arguments:
        * Phase shift gate
        * Rotational gates: (Rx, Ry, Rz)
        * Ising Coupling gates: (Rxx, Ryy, Rzz)
    * Control gates:
        * Ability to easily create a controlled version of an existing gate
        * Multi-controlled and multi-targeted gates
    * Custom gates:
        * Gates created from the current gates in a Sim
        * For black boxes / oracles
        * When printing as a gate within a sim, will print as a single gate
            * A printout similar to the Sim printout can be printed separately
* Support custom qubit ordering for many-bit gates
* Doesn't need to support simulating measurements
    * Reason: Doesn't add much to the project
* Printing features:
    * Print the circuit diagram
        * One qubit per row
        * One gate per column
        * Differentiate control bits from the controlled bit(s)
    * Print current statevector or probability distribution
* Miscellaneous
    * Sims can be named
    * Create a document detailing usage

## Testing

Types of testcases:
    * Gate matrices (base, full) are correct when created with Gates
    * Gate with arguments
    * Controlled gates
    * Custom gates
    * Reordering matrix to LSB order
    * Circuit printouts
    * Statevector consistent with probability distribution in a Sim
    * Correct statevector throughout a Sim

## Future Work

* Support for classical bits and gates
* Support conjugate transpose of a gate matrix as a modifier