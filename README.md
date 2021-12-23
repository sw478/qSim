# Classical Quantum Circuit Simulator and Examples

* Based off the Qiskit API but does not use Qiskit or any other quantum related libraries
* info.md: Details on how to use program
* design_doc: Requirements and types of testcases

## Demos

* Bell State

## Features

* Supports basic quantum gates similar to Qiskit
    * Supports most gates listed from this [Wikipedia](https://en.wikipedia.org/wiki/Quantum_logic_gate) page
* Supports custom qubit ordering for many-bit gates
    * Ex. Toffoli gate with control qubits Q4, Q1 and target qubit Q2
* Printing features:
    * Circuit diagram
    * Current statevector or probability distribution
    * Can inspect (print) between gates
* Doesn't support:
    * Measurements
    * Classical bits