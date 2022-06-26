import numpy as np
from dataclasses import dataclass


@dataclass
class GateContainer:

    """
        gate_label: character representation of gate, printed when circuit diagram is printed

        mat: gate matrix in LSB order, includes control bits

        num_cb: number of control bits for this gate, used for printing

        expected_qbits: inferred expected number of qubits
    """

    gate_label: str
    mat: np.ndarray
    num_cb: int

    def __post_init__(self):
        self.expected_qbits = int(np.log2(self.mat.shape[0]))
