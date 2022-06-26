## quantSim: A Classical Quantum Circuit Simulator with Examples

* Based off the Qiskit API but does not use Qiskit or any other quantum related libraries
* system_design_doc.md: Program usage, architecture, class descriptions, details on error checking
* design_doc.md: Requirements and testcases
* info.md: Usage, list of available gates, miscellaneous project info

## Features

* Supports basic quantum gates similar to Qiskit
    * Supports most gates listed from this [Wikipedia](https://en.wikipedia.org/wiki/Quantum_logic_gate) page
* Supports custom qubit ordering for many-bit gates
    * Ex. Toffoli gate with control qubits Q4, Q1 and target qubit Q2
* Printing features:
    * Circuit diagram
    * Current statevector or probability distribution
    * Can inspect (print) between gate additions
* Doesn't support:
    * Measurements
    * Classical bits
    
## Demos

* Bell State
* Grover's Search Algorithm