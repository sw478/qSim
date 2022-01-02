## Classical Quantum Circuit Simulator and Examples

* Based off the Qiskit API but does not use Qiskit or any other quantum related libraries
* system_design_doc: Program usage, architecture, class descriptions, details on error checking
* design_doc: Requirements and testcases
* info.md: Miscellaneous project info, list of available gates

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
    
## Demos

* Bell State
* Grover's Search Algorithm