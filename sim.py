import numpy as np
from gate import SubGate, GateContainer


class Sim:

    # initializes the simulator with n_qbits qubits
    def __init__(self, n_qbits, name="", for_testing=False):
        self.n_qbits = n_qbits
        self.name = name
        self.for_testing = for_testing
        self.sub_gates = []
        self.as_gate = np.eye(1 << n_qbits, dtype=complex)

    def __repr__(self):
        self.print_sim()

    # adds gate and evaluates the resulting sim as a matrix
    def add_gate(self, gate, i_qbits):
        sub_gate = SubGate(gate, i_qbits, self.n_qbits)
        self.as_gate = np.matmul(sub_gate.full_mat, self.as_gate)

    def to_gate(self, gate_label="G"):
        return GateContainer(gate_label, self.as_gate, 0)

    # prints the quantum algorithm for this sim
    def print_sim(self):
        print(f"[Circuit Diagram] {self.name}")

        for i_qbit in range(self.n_qbits):
            print(f"Q{i_qbit.__repr__()}: ", end="")
            for t, sub_gate in enumerate(self.sub_gates):
                sub_gate.print_gate(i_qbit)

                if t is not len(self.sub_gates) - 1:
                    print("|", end="")
            print("")

        print("")

    def print_statevector(self):
        statevector = np.zeros(1 << self.n_qbits, dtype=complex)
        statevector[0] = 1
        statevector = np.matmul(self.as_gate, statevector)
        statevector = np.around(statevector, 5)

        if not self.for_testing:
            print(f"[Statevector] {self.name}\n{statevector}", end="\n\n")

        return statevector

    def print_prob_dist(self):
        statevector = np.zeros(1 << self.n_qbits, dtype=complex)
        statevector[0] = 1
        statevector = np.matmul(self.as_gate, statevector)
        prob_dist = [np.abs(i) ** 2 for i in statevector]
        prob_dist = np.round(prob_dist, 5)

        if not self.for_testing:
            print(f"[Probability Dist.] {self.name}\n{prob_dist}", end="\n\n")

        return prob_dist
