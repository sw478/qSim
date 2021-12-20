import numpy as np
from gates import Gates, Gate_Inst

"""
    All qubits are initialized to the 0 state [1, 0] at the start

    t: int
    i_qbit: int
    i_qbits: int array
    gate: gate label (GATE_*)
"""


class Sim:
    """
        statevector of qubits

        [i_qbit] = state
    """

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
        returns gate matrix or -1 on error
    """

    def add_gate(self, i_qbits, gate_label, print_mat=False):
        # check values in i_qbits
        if not Gates.check_valid_i_qbits(i_qbits, self.n_qbits):
            return -1

        # check i_qbits is a set
        if not Gates.check_set(i_qbits):
            return -1

        self.gates += [[i_qbits, gate_label]]

        mat = self.__get_mat(i_qbits, gate_label)
        self.qbits = np.matmul(self.qbits, mat)

        if print_mat:
            print(mat)

        return mat

    """
        get the matrix needed for multiplication
    """

    def __get_mat(self, i_qbits, gate_label):
        base_mat = Gates.gate_dict[gate_label]
        gate_size = len(i_qbits)

        if not Gates.check_size(base_mat, i_qbits):
            return -1

        # initial matrix state
        mat = [1]
        for i in range(self.n_qbits - gate_size):
            mat = np.kron(mat, np.eye(2))
        mat = np.kron(mat, base_mat)

        # keep track of i_qbits already swapped
        swap_goal = [i for i in range(self.n_qbits)]
        swap_list = [i for i in range(self.n_qbits)]
        for i in i_qbits:
            swap_list.remove(i)
        swap_list = i_qbits + swap_list

        # loop while list is not in order
        while True:
            # check if list is in order
            for i in range(len(swap_goal)):
                if swap_goal[i] is not swap_list[i]:
                    Gate_Inst.reorder_matrix(mat, i, swap_list[i])

                    # swap list values
                    temp = swap_list[i]
                    swap_list[i] = swap_list[temp]
                    swap_list[temp] = temp
                    continue

            break

        return mat

    """
        prints the quantum algorithm for this sim
    """

    def print_sim(self):
        print("[Circuit Diagram] " + self.name)

        for i_qbit in range(self.n_qbits):
            print("Q" + i_qbit.__repr__() + ": ", end="")
            for t, gate in enumerate(self.gates):
                i_qbits = gate[0]
                gate_label = gate[1]
                Gates.print_gate(gate_label, i_qbits, i_qbit)

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
