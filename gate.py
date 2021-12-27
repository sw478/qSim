import numpy as np
from dataclasses import dataclass, field

from typing import List

"""
    Instance of a gate in a Sim or a custom gate
    i_qbits: mapping of Sim's qubit indices to gate's qubit indices
    n_qbits: number of qubits in Sim, determines size of full gate matrix
    full matrix: tensor product of gate matrix and identity gates in parallel
"""


class SubGate:
    full_mat = None

    def __init__(self, gate, i_qbits, n_qbits):
        self.gate = gate
        self.i_qbits = i_qbits
        self.n_qbits = n_qbits

        if not Error.check_size(gate.mat, i_qbits):
            return

        # check values in i_qbits
        if not Error.check_valid_i_qbits(i_qbits, self.n_qbits):
            return

        # check i_qbits is a set
        if not Error.check_set(i_qbits):
            return

        self.full_mat = self.__get_full_mat()

    def __get_full_mat(self):
        # initial unswapped full matrix in LSB
        mat = [1]
        for i in range(self.n_qbits - len(self.i_qbits)):
            mat = np.kron(mat, np.eye(2))
        mat = np.kron(mat, self.gate.mat)

        # goal permutation
        swap_goal = [i for i in range(self.n_qbits)]

        # current permutation
        gen = (i for i in swap_goal if i not in self.i_qbits)
        swap_list = self.i_qbits + [i for i in gen]

        # permute until goal is reached
        for i, (a, b) in enumerate(zip(swap_list, swap_goal)):
            if b == a:
                continue

            full_swap_mat = Permute.get_swap_mat(self.n_qbits, i, a)
            mat = np.linalg.multi_dot([full_swap_mat, mat, full_swap_mat])
            swap_list[i], swap_list[a] = swap_list[a], swap_list[i]

        return mat

    def print_gate(self, i_qbit):
        num_chars = len(self.gate.gate_label)

        if i_qbit not in self.i_qbits:
            [print(" ", end="") for i in range(num_chars)]
            if self.gate.num_cb > 0:
                print(" ", end="")
            return

        j_qbit = self.i_qbits.index(i_qbit)
        if j_qbit < self.gate.num_cb:
            print("C", end="")
            [print(" ", end="") for i in range(num_chars)]
        else:
            if self.gate.num_cb > 0:
                print(" ", end="")
            print(self.gate.gate_label, end="")


"""
    gate_label: character representation of gate, printed when circuit diagram is printed
    mat: gate matrix in LSB order, includes control bits
    num_cb: number of control bits for this gate, only used for printing
    sub_gates: if this is a custom gate, then this will be the list of SubGates in a Sim
"""


@dataclass
class GateContainer:
    gate_label: str
    mat: np.ndarray
    num_cb: int

    def __post_init__(self):
        self.expected_qbits = int(np.log2(self.mat.shape[0]))


class Permute:
    swap_mat = np.array(
        [[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1]])

    __max_n_qbits = 0
    __preproc = None

    @staticmethod
    def get_swap_mat(n_qbits, i_qbit, j_qbit):
        if Permute.__max_n_qbits < n_qbits:
            Permute.preprocess(n_qbits)

        if i_qbit > j_qbit:
            i_qbit, j_qbit = j_qbit, i_qbit

        full_swap_mat = Permute.__preproc[i_qbit][j_qbit]

        if Permute.__max_n_qbits > n_qbits:
            mat_size = 1 << n_qbits
            full_swap_mat = full_swap_mat[0:mat_size, 0:mat_size]

        return full_swap_mat

    @staticmethod
    def preprocess(n_qbits):
        mat_size = 1 << n_qbits
        Permute.__preproc = np.ndarray(shape=(n_qbits, n_qbits, mat_size, mat_size))
        Permute.__max_n_qbits = n_qbits

        # [i][i+1] diagonal
        for i in range(n_qbits - 1):
            mat = [1]

            for k in range(n_qbits):
                if k == i + 1:
                    continue
                elif k == i:
                    mat = np.kron(Permute.swap_mat, mat)
                else:
                    mat = np.kron(np.eye(2), mat)

            Permute.__preproc[i][i + 1] = mat
            Permute.__preproc[i + 1][i] = mat

        for j in range(2, n_qbits):
            for i in range(j - 2, -1, -1):
                a = Permute.__preproc[i][i + 1]
                b = Permute.__preproc[i + 1][j]
                Permute.__preproc[i][j] = np.linalg.multi_dot([a, b, a])


"""
    Gate definitions
"""


class Gate:
    SQRT_2 = 1 / np.sqrt(2)
    SQRT_8 = 1 / np.sqrt(8)
    SQRT_32 = 1 / np.sqrt(32)

    @staticmethod
    def I():
        mat = np.eye(2)
        return GateContainer("I", mat, 0)

    @staticmethod
    def H():
        mat = np.array(
            [[Gate.SQRT_2, Gate.SQRT_2],
             [Gate.SQRT_2, -Gate.SQRT_2]])

        return GateContainer("H", mat, 0)

    @staticmethod
    def X():
        mat = np.array(
            [[0, 1],
             [1, 0]])

        return GateContainer("X", mat, 0)

    @staticmethod
    def Y():
        mat = np.array(
            [[0, -1j],
             [1j, 0]])

        return GateContainer("Y", mat, 0)

    @staticmethod
    def Z():
        mat = np.array(
            [[1, 0],
             [0, -1]])

        return GateContainer("Z", mat, 0)

    # sqrt X
    @staticmethod
    def rX():
        mat = np.array(
            [[0.5 + 0.5j, 0.5 - 0.5j],
             [0.5 - 0.5j, 0.5 + 0.5j]])

        return GateContainer("rX", mat, 0)

    @staticmethod
    def P(theta, gate_label="P"):
        mat = np.array(
            [[1, 0],
             [0, np.exp(theta * 1j)]])

        return GateContainer(gate_label, mat, 0)

    @staticmethod
    def Rx(theta):
        theta_2 = theta * 0.5
        mat = np.array(
            [[np.cos(theta_2), np.sin(theta_2) * -1j],
             [np.sin(theta_2) * -1j, np.cos(theta_2)]])

        return GateContainer("Rx", mat, 0)

    @staticmethod
    def Ry(theta):
        theta_2 = theta * 0.5
        mat = np.array(
            [[np.cos(theta_2), np.sin(theta_2) * -1],
             [np.sin(theta_2), np.cos(theta_2)]])

        return GateContainer("Ry", mat, 0)

    @staticmethod
    def Rz(theta):
        theta_2 = theta * 0.5
        mat = np.array(
            [[np.exp(-1j * theta_2), 0],
             [0, np.exp(1j * theta_2)]])

        return GateContainer("Rz", mat, 0)

    @staticmethod
    def Rxx(theta):
        theta_2 = theta * 0.5
        cos_theta_2 = np.cos(theta_2)
        isin_theta_2_neg = -1j * np.sin(theta_2)

        mat = np.array(
            [[cos_theta_2, 0, 0, isin_theta_2_neg],
             [0, cos_theta_2, isin_theta_2_neg, 0],
             [0, isin_theta_2_neg, cos_theta_2, 0],
             [isin_theta_2_neg, 0, 0, cos_theta_2]])

        return GateContainer("Rxx", mat, 0)

    @staticmethod
    def Ryy(theta):
        theta_2 = theta * 0.5
        cos_theta_2 = np.cos(theta_2)
        isin_theta_2_neg = -1j * np.sin(theta_2)

        mat = np.array(
            [[cos_theta_2, 0, 0, -isin_theta_2_neg],
             [0, cos_theta_2, isin_theta_2_neg, 0],
             [0, isin_theta_2_neg, cos_theta_2, 0],
             [-isin_theta_2_neg, 0, 0, cos_theta_2]])

        return GateContainer("Ryy", mat, 0)

    @staticmethod
    def Rzz(theta):
        theta_2 = theta * 0.5
        e_theta_2 = np.exp(1j * theta_2)
        e_theta_2_neg = np.exp(-1j * theta_2)

        mat = np.array(
            [[e_theta_2_neg, 0, 0, 0],
             [0, e_theta_2, 0, 0],
             [0, 0, e_theta_2, 0],
             [0, 0, 0, e_theta_2_neg]])

        return GateContainer("Rzz", mat, 0)

    """
        Gates created with other matrices
    """

    @staticmethod
    def S(dagger=False):
        if dagger:
            return Gate.Sdg()

        mat = np.array(
            [[1, 0],
             [0, 1j]])

        return GateContainer("S", mat, 0)

    @staticmethod
    def Sdg():
        mat = np.array(
            [[1, 0],
             [0, -1j]])

        return GateContainer("Sdg", mat, 0)

    @staticmethod
    def T(dagger=False):
        if dagger:
            return Gate.Tdg()

        mat = np.array(
            [[1, 0],
             [0, Gate.SQRT_2 + 1j * Gate.SQRT_2]])

        return GateContainer("T", mat, 0)

    @staticmethod
    def Tdg():

        mat = np.array(
            [[1, 0],
             [0, Gate.SQRT_2 - 1j * Gate.SQRT_2]])

        return GateContainer("Tdg", mat, 0)

    @staticmethod
    def CX():
        gate = Gate.X()
        return Gate.MC(gate, num_cb=1)

    @staticmethod
    def Toffoli():
        gate = Gate.X()
        return Gate.MC(gate, num_cb=2)

    @staticmethod
    def CSwap():
        gate = Gate.Swap()
        return Gate.MC(gate, num_cb=1)

    @staticmethod
    def CP(theta):
        gate = Gate.P(theta)
        return Gate.MC(gate, num_cb=1)

    @staticmethod
    def Swap():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [0, 0, 0, 1]])

        return GateContainer("SW", mat, 0)

    @staticmethod
    def Swap_sqrt():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0.5 + 0.5j, 0.5 - 0.5j, 0],
             [0, 0.5 - 0.5j, 0.5 + 0.5j, 0],
             [0, 0, 0, 1]])

        return GateContainer("rSW", mat, 0)

    @staticmethod
    def Swap_i():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1j, 0],
             [0, 1j, 0, 0],
             [0, 0, 0, 1]])

        return GateContainer("iSW", mat, 0)

    @staticmethod
    def Swap_sqrt_i():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, Gate.SQRT_2, 1j * Gate.SQRT_2, 0],
             [0, 1j * Gate.SQRT_2, Gate.SQRT_2, 0],
             [0, 0, 0, 1]])

        return GateContainer("irSW", mat, 0)

    """
        Special gate operations
    """

    # conjugate transpose of matrix
    @staticmethod
    def dg(gate):
        mat = np.matrix.h(gate.mat)

        return GateContainer(gate.gate_label, mat, gate.num_cb)

    # multi-targeted gate for single bit gates
    # num_tb: number of target bits
    @staticmethod
    def MT(gate, num_tb):

        n_qbits_mat = np.log2(gate.mat.shape[0])
        if n_qbits_mat != 1:
            print(f"Error: Gate is not single bit: {n_qbits_mat}")
            return

        mat = [1]
        for i in range(num_tb):
            mat = np.kron(mat, gate.mat)

        return GateContainer(gate.gate_label, mat, 0)

    # multi-controlled gate
    # num_cb: number of control bits
    # recursively creates a controlled version of gate.mat
    @staticmethod
    def MC(gate, num_cb):
        control_mat = np.array([[1, 0], [0, 0]])
        target_mat = np.array([[0, 0], [0, 1]])

        mat = np.copy(gate.mat)
        eye_size = mat.shape[0]
        for i in range(num_cb):
            mat = np.kron(mat, target_mat)
            eye = np.eye(eye_size)
            eye_size *= 2
            mat = mat + np.kron(eye, control_mat)

        return GateContainer(gate.gate_label, mat, gate.num_cb + num_cb)


class Error:

    @staticmethod
    def check_size(mat, i_qbits):
        len_i_qbits = len(i_qbits)
        len_matrix = np.log2(mat.shape[0])

        if len_i_qbits != len_matrix:
            print(f"Error: Dimension mismatch: {str(len_i_qbits)} {str(len_matrix)}")
            return False
        return True

    @staticmethod
    def check_set(i_qbits):
        if len(i_qbits) is not len(list(dict.fromkeys(i_qbits))):
            print(f"Error: i_qbits is not a set: {i_qbits.__repr__()}")
            return False

        return True

    @staticmethod
    def check_valid_i_qbits(i_qbits, n_qbits):
        for i_qbit in i_qbits:
            if i_qbit >= n_qbits:
                print(f"Error: Invalid i_qbits values: {i_qbits.__repr__()}")
                return False

        return True
