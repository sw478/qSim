
### Project Info

* LSB bit ordering similar to [Qiskit](https://qiskit.org/documentation/tutorials/circuits/3_summary_of_quantum_operations.html#Basis-vector-ordering-in-Qiskit)
* Gates with control bits have the control bits as the MSBs
* All qubits are initialized to [0, 1].

#### Currently Implemented Gates

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