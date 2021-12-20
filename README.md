# Classical Quantum Circuit Simulator and Examples

## General Notes

* Based off the Qiskit API but does not use Qiskit or any quantum libaries
* LSB bit ordering similar to [Qiskit](https://qiskit.org/documentation/tutorials/circuits/3_summary_of_quantum_operations.html#Basis-vector-ordering-in-Qiskit)
* Gate matrices are multiplied once gates are added, so the statevector can be viewed in between gate additions.
    I am aware this [isn't possible](https://en.wikipedia.org/wiki/Measurement_in_quantum_mechanics)
    in an actual quantum computer, but it is designed this way for debugging/educational purposes.

## Functionalities

* Supports basic quantum gates
* Supports custom qubit ordering for many-bit gates
    * Example: can support a toffoli gate with control qubits at 4, 1 and controlled qubit at 0
* Doesn't support the idea of measurements, but the statevector and probability distributions can be printed
* Doesn't support classical bits
* Displays gates in a circuit diagram, can show which qubits are control/controlled

## Simulation Examples

* Bell State
* CCC-Not universal simplification with Toffoli gates
