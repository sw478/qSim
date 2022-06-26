import numpy as np
from quantSim import Error
from quantSim import Permute


class SubGate:

    full_mat = None

    def __init__(self, gate, i_qbits, n_qbits):
        """
        Instance of a gate in a Sim or a custom gate

        full_mat: tensor product of gate matrix and identity gates in parallel

        :param gate: gate with potentially fewer qubits than Sim has
        :param i_qbits: mapping of Sim's qubit indices to gate's qubit indices
        :param n_qbits: number of qubits in Sim, determines size of full gate matrix
        """
        self.gate = gate
        self.__i_qbits = i_qbits
        self.__n_qbits = n_qbits

        if not Error.check_size(gate.mat, i_qbits):
            return

        if not Error.check_valid_i_qbits(i_qbits, self.__n_qbits):
            return

        if not Error.check_set(i_qbits):
            return

        self.full_mat = self.__get_full_mat()

    def __get_full_mat(self):
        """
        Returns the gate represented within the Sim
        """
        # initial unswapped full matrix in LSB
        mat = [1]
        for i in range(self.__n_qbits - len(self.__i_qbits)):
            mat = np.kron(mat, np.eye(2))
        mat = np.kron(mat, self.gate.mat)

        # goal permutation
        swap_goal = [i for i in range(self.__n_qbits)]

        # current permutation
        gen = (i for i in swap_goal if i not in self.__i_qbits)
        swap_list = self.__i_qbits + [i for i in gen]

        # permute until goal is reached
        for i, (a, b) in enumerate(zip(swap_list, swap_goal)):
            if b == a:
                continue

            full_swap_mat = Permute.get_swap_mat(self.__n_qbits, i, a)
            mat = np.linalg.multi_dot([full_swap_mat, mat, full_swap_mat])
            swap_list[i], swap_list[a] = swap_list[a], swap_list[i]

        return mat

    # i_qbit referring to the
    def print_gate(self, i_qbit):
        """
        Prints gate label for sim visualization

        :param i_qbit: Sim i_qbit
        :return: None
        """
        num_chars = len(self.gate.gate_label)

        if i_qbit not in self.__i_qbits:
            [print(" ", end="") for i in range(num_chars)]
            if self.gate.num_cb > 0:
                print(" ", end="")
            return

        j_qbit = self.__i_qbits.index(i_qbit)
        if j_qbit < self.gate.num_cb:
            print("C", end="")
            [print(" ", end="") for i in range(num_chars)]
        else:
            if self.gate.num_cb > 0:
                print(" ", end="")
            print(self.gate.gate_label, end="")
