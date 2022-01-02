## Design Doc

### Problem Description

The goal of this project is to create a classical quantum circuit simulator.
This will be mainly used as practice for coding and as a learning experience for quantum computing.

### Requirements

* Create quantum circuits similar to online quantum computer simulators
* Gates to support:
    * Basic gates:
        * Identity
        * Pauli (X, Y, Z)
        * Hadamard
        * Phase shift: (S, T, Sdg, Tdg)
        * Swap and variations: (sqrt, imaginary)
    * Multi-targeted single bit gates
    * Multi-controlled gates
    * Gates with arguments:
        * Phase shift gate
        * Rotational gates: (Rx, Ry, Rz)
        * Ising Coupling gates: (Rxx, Ryy, Rzz)
    * Custom gates:
        * Single gate created from the current gates in a Sim
* Support custom qubit ordering for many-bit gates
* Doesn't need to support simulating measurements
    * Reason: Doesn't add much to the project
* Printing features:
    * Print the circuit diagram
        * One qubit per row
        * One gate per column
        * Differentiate control bits from the target bit(s)
        * When printing a custom gate within a printout, print as a single gate
    * Print current statevector or probability distribution

### Testing

Types of testcases:
* Gate matrices (base, full) are correct when created with Gates
* Gate with arguments
* Control gates
* Custom gates
* Reordering matrix to LSB order
* Circuit printouts
* Statevector consistent with probability distribution in a Sim
* Correct statevector throughout a Sim

### Error checking

List of common things needed to be checked:
* check_size(): ensures the size of a matrix corresponds to the correct number of qubits
* check_set(): ensures "i_qbits" lists are a set
* check_valid_i_qbits(): ensures "i_qbits" values are valid (lower than n_qbits)

### Future Work

* Support for classical bits and gates
* Support custom initial statevector in Sims
* Enable reuse of array of swap matrices when growing
* Add more gates and demos