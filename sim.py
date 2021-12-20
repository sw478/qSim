import numpy as np
from gates import Gates, Gate_Inst

"""
    All qubits are initialized to the 0 state [1, 0] at the start

    i_qbit: int
    i_qbits: int array
    gate: gate label (GATE_*)
"""


class Sim:

    # statevector of qubits
    qbits = np.array([], dtype=complex)

    """
        initializes the simulator with n_qbits qubits
            and n_cbits classical bits
    """

    def __init__(self, n_qbits, name="", for_testing=False):
        self.n_qbits = n_qbits
        self.name = name
        self.for_testing = for_testing

        """
            list of gates
            entries: [i_qbits, gate]
        """
        self.gates = []

        """
            statevector of all qubits
            initialize all qubits to 0
        """
        self.qbits = np.zeros(1 << n_qbits, dtype=complex)
        self.qbits[0] = 1

    """
        adds the gate to the simulator with the qbits listed in i_qbits
        evaluates the resulting matrix with the statevector
    """

    def add_gate(self, gate_inst):
        self.gates += [gate_inst]
        self.qbits = np.matmul(self.qbits, gate_inst.full_mat)

    """
        prints the quantum algorithm for this sim
    """

    def print_sim(self):
        print("[Circuit Diagram] " + self.name)

        for i_qbit in range(self.n_qbits):
            print("Q" + i_qbit.__repr__() + ": ", end="")
            for t, gate_inst in enumerate(self.gates):
                Gates.print_gate(gate_inst, i_qbit)

                if t is not len(self.gates) - 1:
                    print("|", end="")
            print("")

        print("")

    def print_statevector(self):
        if not self.for_testing:
            print("[Statevector] " + self.name)
            print(self.qbits, end="\n\n")

        return self.qbits

    def print_prob_dist(self):
        prob_dist = (self.qbits ** 2).real

        if not self.for_testing:
            print("[Probability Dist.] " + self.name)
            print(prob_dist, end="\n\n")

        return prob_dist
