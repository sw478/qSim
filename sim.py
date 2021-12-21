import numpy as np
from gates import Gates, Gate_Inst, Gate

"""
    All qubits are initialized to the 0 state [1, 0] at the start

    i_qbit: int
    i_qbits: int array
    gate: gate label (GATE_*)
"""


class Sim:

    # initializes the simulator with n_qbits qubits
    def __init__(self, n_qbits, name="", for_testing=False):
        self.n_qbits = n_qbits
        self.name = name
        self.for_testing = for_testing

        # list of gates
        self.gates = []

        # matrix if turning sim into gate
        self.as_gate = np.eye(1 << n_qbits, dtype=complex)

    # adds the gate and evaluates the resulting sim as a matrix
    def add_gate(self, gate, i_qbits):
        gate_inst = Gate_Inst(gate, i_qbits, self.n_qbits)
        self.gates += [gate_inst]
        self.as_gate = np.matmul(self.as_gate, gate_inst.full_mat)

    def turn_into_gate(self, gate_label="G"):
        return Gate(gate_label, self.as_gate, 0, self.gates)

    # prints the quantum algorithm for this sim
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
        statevector = np.zeros(1 << self.n_qbits, dtype=complex)
        statevector[0] = 1
        statevector = np.matmul(statevector, self.as_gate)

        if not self.for_testing:
            print("[Statevector] " + self.name)
            print(statevector, end="\n\n")

        return statevector

    def print_prob_dist(self):
        statevector = np.zeros(1 << self.n_qbits, dtype=complex)
        statevector[0] = 1
        statevector = np.matmul(statevector, self.as_gate)
        prob_dist = (statevector ** 2).real

        if not self.for_testing:
            print("[Probability Dist.] " + self.name)
            print(prob_dist, end="\n\n")

        return prob_dist
