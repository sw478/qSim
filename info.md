## Project Info

* LSB bit ordering similar to [Qiskit](https://qiskit.org/documentation/tutorials/circuits/3_summary_of_quantum_operations.html#Basis-vector-ordering-in-Qiskit)
* Gates with control bits have the control bits as the MSBs
* All qubits are initialized to [0, 1].

### Currently Implemented Gates

* Identity
* Pauli (X, Y, Z)
* Hadamard
* Sqrt X (rX)
* Phase shift: (S, T, Sdg, Tdg)
* Swap and variations: (sqrt, imaginary)
* Multi-targeted single bit gates
* Multi-controlled gates
* Gates with arguments:
    * Phase shift gate
    * Rotational gates: (Rx, Ry, Rz)
    * Ising Coupling gates: (Rxx, Ryy, Rzz)
* Custom gates:
    * Single black box gate created from the current gates in a Sim

### Usage

#### Creating a simulation

1. Instantiate a Sim with the desired number of qubits.
2. Add gates sequentially with add_gate(gate, i_qbits).
    a. gate: Gate object describing the gate used
    b. i_qbits: Ordered set of indices of qubits to be used for the gate
3. Call print_statevector(), print_prob_dist(), print_sim() to print Sim info.

#### Creating a custom gate from a Sim

1. Create a Sim with desired gates
2. Call to_gate() to return a Gate object of the Sim
3. Can now use Gate in another Sim