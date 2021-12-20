import numpy as np
from gates import Gates

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
        list of gates (np arrays) of type "GATE_*"

        [t] = [i_qbits, gate]
    """

    gates = {}

    """
        initializes the simulator with n_qbits qubits
            and n_cbits classical bits
    """

    def __init__(self, n_qbits, name="", for_testing=False):
        self.n_qbits = n_qbits
        self.space_size = 1 << n_qbits
        self.name = name
        self.for_testing = for_testing
        self.max_t = 0

        # initialize all qubits to |0>
        self.qbits = np.zeros(2 ** self.n_qbits, dtype=complex)
        self.qbits[0] = 1

    """
        adds the gate to the simulator with the qbits listed in i_qbits
        returns gate matrix or -1 on error
    """

    def add_gate(self, i_qbits, gate_label, print_mat=False):
        # check values in i_qbits
        if self.__check_valid_i_qbits(i_qbits) != 0:
            print("Cannot add gate, invalid i_qbits values: " + i_qbits.__repr__())
            return -1

        # check i_qbits is a set
        if len(i_qbits) is not len(list(dict.fromkeys(i_qbits))):
            print("Cannot add gate, i_qbits is not a set: " + i_qbits.__repr__())
            return -1

        self.gates[self.max_t] = [i_qbits, gate_label]
        self.max_t += 1

        mat = self.__get_mat(i_qbits, gate_label)
        self.qbits = np.matmul(self.qbits, mat)

        if print_mat:
            print(mat)

        return mat

    def __check_valid_i_qbits(self, i_qbits):
        for i_qbit in i_qbits:
            if i_qbit >= self.n_qbits:
                return -1

        return 0

    """
        get the matrix needed for multiplication
    """

    def __get_mat(self, i_qbits, gate_label):
        base_mat = Gates.gate_dict[gate_label]
        gate_size = len(i_qbits)

        if gate_size != np.log2(base_mat.shape[0]):
            print("dimension mismatch: i_qbits: " + i_qbits.__repr__())
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
                    self.reorder_matrix(mat, i, swap_list[i])

                    # swap list values
                    temp = swap_list[i]
                    swap_list[i] = swap_list[temp]
                    swap_list[temp] = temp
                    continue

            break

        return mat

    """
        swaps rows/cols of matrix to match gate's i_qbits
        enables custom qubit ordering for many-bit gates
    """

    def reorder_matrix(self, mat, i_qbit, j_qbit):
        if i_qbit == j_qbit:
            return

        #print("swapping qubits: [" + str(i_qbit) + " " + str(j_qbit) + "]")

        # keep track of rows/cols already swapped
        swap_list = [i for i in range(self.space_size)]

        # ensures i < j for positive bit shift value
        if i_qbit > j_qbit:
            i_qbit, j_qbit = j_qbit, i_qbit
        shift = j_qbit - i_qbit

        # get row/col index (bit)
        i2 = 1 << i_qbit
        j2 = 1 << j_qbit

        while swap_list:
            # a & b are matrix row/col indices
            a = swap_list[0]
            a_i = a & i2
            a_j = a & j2

            if a_i is (a_j >> shift):
                swap_list.remove(a)
                continue

            # b: a if bits at i2 and j2 were swapped
            b = a & ~(i2 | j2)  # set bits i2 and j2 to 0
            b |= a_j >> shift  # set bit i2 to bit value j2
            b |= a_i << shift  # set bit j2 to bit value i2

            # check if matrix needs to be reordered
            if a ^ b:
                # print("[" + str(i_qbit) + " " + str(j_qbit) + "] swapping: " + str(a) + " " + str(b))
                mat[[a, b]] = mat[[b, a]]  # swap rows
                mat[:, [a, b]] = mat[:, [b, a]]  # swap columns

            swap_list.remove(a)
            swap_list.remove(b)

    """
        prints the quantum algorithm for this sim
    """

    def print_sim(self):
        print("[Circuit Diagram] " + self.name)

        for i_qbit in range(self.n_qbits):
            print("Q" + i_qbit.__repr__() + ": ", end="")
            for t in range(self.max_t):
                gate = self.gates[t]

                i_qbits = gate[0]
                gate_label = gate[1]
                Gates.print_gate(gate_label, i_qbits, i_qbit)

                if t is not self.max_t - 1:
                    print("|", end="")
            print("")

        print("")

    def print_statevector(self):

        if self.for_testing:
            return self.qbits

        print("[Statevector] " + self.name)
        print(self.qbits, end="\n\n")

    def print_prob_dist(self):
        prob_dist = (self.qbits ** 2).real

        if self.for_testing:
            return prob_dist

        print("[Probability Dist.] " + self.name)
        print(prob_dist, end="\n\n")
