import numpy as np

"""
    gate_label: characters to be printed when circuit diagram is printed (not including the 'C's for control bits)
    mat: gate matrix
        -LSB
        -includes control bits
    num_cb: number of control bits for this gate, control bits are the MSBs, only used for printing
"""


class Gate:

    def __init__(self, gate_label, mat, num_cb):
        self.gate_label = gate_label
        self.mat = mat
        self.num_cb = num_cb

"""
    Instance of a gate in a Sim or a custom gate
    gate: Gate object
    i_qbits: mapping of Sim's qubit indices to gate's qubit indices
    n_qbits: number of qubits in Sim, determines size of full gate matrix
"""


class Gate_Inst:
    full_mat = None

    def __init__(self, gate, i_qbits, n_qbits):
        self.gate = gate
        self.i_qbits = i_qbits
        self.n_qbits = n_qbits

        if not Gates.check_size(gate.mat, i_qbits):
            return

        # check values in i_qbits
        if not Gates.check_valid_i_qbits(i_qbits, self.n_qbits):
            return

        # check i_qbits is a set
        if not Gates.check_set(i_qbits):
            return

        self.full_mat = self.__get_full_mat()

    def __get_full_mat(self):
        # initial unswapped full matrix
        mat = [1]
        for i in range(self.n_qbits - len(self.i_qbits)):
            mat = np.kron(mat, np.eye(2))
        mat = np.kron(mat, self.gate.mat)

        # keep track of i_qbits already swapped
        swap_goal = [i for i in range(self.n_qbits)]
        swap_list = [i for i in range(self.n_qbits)]
        for i in self.i_qbits:
            swap_list.remove(i)
        swap_list = self.i_qbits + swap_list

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
        swaps rows/cols of matrix to match gate's i_qbits
        enables custom qubit ordering for many-bit gates
    """

    @staticmethod
    def reorder_matrix(mat, i_qbit, j_qbit):
        if i_qbit == j_qbit:
            return

        n_elements = mat.shape[0]

        if i_qbit | j_qbit > np.log2(n_elements):
            return

        #print("swapping qubits: [" + str(i_qbit) + " " + str(j_qbit) + "]")

        # keep track of rows/cols already swapped
        swap_list = [i for i in range(n_elements)]

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


class Gates:
    SQRT_H = 1 / np.sqrt(2)
    SQRT_E = 1 / np.sqrt(8)

    @staticmethod
    def I(i_qbits, n_qbits):
        gate = Gates.I_mat()
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def H(i_qbits, n_qbits):
        gate = Gates.H_mat()
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def X(i_qbits, n_qbits):
        gate = Gates.X_mat()
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Y(i_qbits, n_qbits):
        gate = Gates.Y_mat()
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Z(i_qbits, n_qbits):
        gate = Gates.Z_mat()
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def P(i_qbits, n_qbits, theta):
        gate = Gates.P_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Rx(i_qbits, n_qbits, theta):
        gate = Gates.Rx_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Ry(i_qbits, n_qbits, theta):
        gate = Gates.Ry_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Rz(i_qbits, n_qbits, theta):
        gate = Gates.Rz_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Rxx(i_qbits, n_qbits, theta):
        gate = Gates.Rxx_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Ryy(i_qbits, n_qbits, theta):
        gate = Gates.Ry_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Rzz(i_qbits, n_qbits, theta):
        gate = Gates.Rz_mat(theta)
        return Gate_Inst(gate, i_qbits, n_qbits)

    """
        Gates created with other matrices
    """
    @staticmethod
    def S(i_qbits, n_qbits, dagger=False):
        theta = -np.pi / 2

        if not dagger:
            gate = Gates.P_mat(theta, "S")
        else:
            gate = Gates.P_mat(-theta, "Sdg")

        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def T(i_qbits, n_qbits, dagger=False):
        theta = -np.pi / 4

        if not dagger:
            gate = Gates.P_mat(theta, "T")
        else:
            gate = Gates.P_mat(-theta, "Tdg")

        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def CX(i_qbits, n_qbits):
        gate = Gates.X_mat()
        gate = Gates.C_mat(gate, num_cb=1)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Toffoli(i_qbits, n_qbits):
        gate = Gates.X_mat()
        gate = Gates.C_mat(gate, num_cb=2)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def CP(i_qbits, n_qbits, theta):
        gate = Gates.P_mat(theta)
        gate = Gates.C_mat(gate, num_cb=1)
        return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def Swap(i_qbits, n_qbits, imag=False, sqrt=False):
        if imag and sqrt:
            gate = Gates.irSwap_mat()
        elif not imag and sqrt:
            gate = Gates.rSwap_mat()
        elif imag and not sqrt:
            gate = Gates.iSwap_mat()
        else:
            gate = Gates.Swap_mat()

        return Gate_Inst(gate, i_qbits, n_qbits)

    # sqrt X
    @staticmethod
    def rX(i_qbits, n_qbits):
        gate = Gates.rX_mat()
        return Gate_Inst(gate, i_qbits, n_qbits)

    """
        Gates generated from a Gate_Inst
    """

    # multi-controlled gate
    @staticmethod
    def MC(gate_inst, ctrl_bits):
        num_cb = len(ctrl_bits)

        if not Gates.check_set(ctrl_bits):
            return

        if not Gates.check_valid_i_qbits(ctrl_bits, gate_inst.n_qbits):
            return

        if num_cb < 1:
            print("Error: Invalid number of control bits:")
            return

        gate = Gates.C_mat(gate_inst.full_mat, num_cb)

        i_qbits = ctrl_bits + gate_inst.i_qbits
        return Gate_Inst(gate, i_qbits, gate_inst.n_qbits)


    # matrix definitions
    @staticmethod
    def I_mat():
        mat = np.eye(2)
        return Gate("I", mat, 0)

    @staticmethod
    def H_mat():
        mat = np.array(
            [[Gates.SQRT_H, Gates.SQRT_H],
             [Gates.SQRT_H, -Gates.SQRT_H]])

        return Gate("H", mat, 0)

    @staticmethod
    def X_mat():
        mat = np.array(
            [[0, 1],
             [1, 0]])

        return Gate("X", mat, 0)

    @staticmethod
    def Y_mat():
        mat = np.array(
            [[0, -1j],
             [1j, 0]])

        return Gate("Y", mat, 0)

    @staticmethod
    def Z_mat():
        return Gates.P_mat(np.pi, "Z")

    @staticmethod
    def S_mat():
        return Gates.P_mat(np.pi / 2, "S")

    @staticmethod
    def Sdg_mat():
        return Gates.P_mat(-np.pi / 2, "Sdg")

    @staticmethod
    def T_mat():
        return Gates.P_mat(np.pi / 4, "T")

    @staticmethod
    def Tdg_mat():
        return Gates.P_mat(np.pi / 4, "Tdg")

    @staticmethod
    def P_mat(theta, gate_label="P"):
        mat = np.array(
            [[1, 0],
             [0, np.exp(theta * 1j)]])

        return Gate(gate_label, mat, 0)

    @staticmethod
    def rX_mat():
        mat = np.array(
            [[0.5 + 0.5j, 0.5 - 0.5j],
             [0.5 - 0.5j, 0.5 + 0.5j]])

        return Gate("rX", mat, 0)

    @staticmethod
    def Swap_mat():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [0, 0, 0, 1]])

        return Gate("SW", mat, 0)

    @staticmethod
    def rSwap_mat():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0.5 + 0.5j, 0.5 - 0.5j, 0],
             [0, 0.5 - 0.5j, 0.5 + 0.5j, 0],
             [0, 0, 0, 1]])

        return Gate("rSW", mat, 0)

    @staticmethod
    def iSwap_mat():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1j, 0],
             [0, 1j, 0, 0],
             [0, 0, 0, 1]])

        return Gate("iSW", mat, 0)

    @staticmethod
    def irSwap_mat():
        mat = np.array(
            [[1, 0, 0, 0],
             [0, Gates.SQRT_H, 1j * Gates.SQRT_H, 0],
             [0, 1j * Gates.SQRT_H, Gates.SQRT_H, 0],
             [0, 0, 0, 1]])

        return Gate("irSW", mat, 0)

    @staticmethod
    def Rx_mat(theta):
        theta_2 = theta * 0.5
        mat = np.array(
            [[np.cos(theta_2), np.sin(theta_2) * -1j],
             [np.sin(theta_2) * -1j, np.cos(theta_2)]])

        return Gate("Rx", mat, 0)

    @staticmethod
    def Ry_mat(theta):
        theta_2 = theta * 0.5
        mat = np.array(
            [[np.cos(theta_2), np.sin(theta_2) * -1],
             [np.sin(theta_2), np.cos(theta_2)]])

        return Gate("Ry", mat, 0)

    @staticmethod
    def Rz_mat(theta):
        theta_2 = theta * 0.5
        mat = np.array(
            [[np.exp(-1j * theta_2), 0],
             [0, np.exp(1j * theta_2)]])

        return Gate("Rz", mat, 0)

    @staticmethod
    def Rxx_mat(theta):
        theta_2 = theta * 0.5
        cos_theta_2 = np.cos(theta_2)
        isin_theta_2_neg = -1j * np.sin(theta_2)

        mat = np.array(
            [[cos_theta_2, 0, 0, isin_theta_2_neg],
             [0, cos_theta_2, isin_theta_2_neg, 0],
             [0, isin_theta_2_neg, cos_theta_2, 0],
             [isin_theta_2_neg, 0, 0, cos_theta_2]])

        return Gate("Rxx", mat, 0)

    @staticmethod
    def Ryy_mat(theta):
        theta_2 = theta * 0.5
        cos_theta_2 = np.cos(theta_2)
        isin_theta_2_neg = -1j * np.sin(theta_2)

        mat = np.array(
            [[cos_theta_2, 0, 0, -isin_theta_2_neg],
             [0, cos_theta_2, isin_theta_2_neg, 0],
             [0, isin_theta_2_neg, cos_theta_2, 0],
             [-isin_theta_2_neg, 0, 0, cos_theta_2]])

        return Gate("Ryy", mat, 0)

    @staticmethod
    def Rzz_mat(theta):
        theta_2 = theta * 0.5
        e_theta_2 = np.exp(1j * theta_2)
        e_theta_2_neg = np.exp(-1j * theta_2)

        mat = np.array(
            [[e_theta_2_neg, 0, 0, 0],
             [0, e_theta_2, 0, 0],
             [0, 0, e_theta_2, 0],
             [0, 0, 0, e_theta_2_neg]])

        return Gate("Rzz", mat, 0)

    # recursively creates a controlled version of gate.mat with num_cb number of control bits
    @staticmethod
    def C_mat(gate, num_cb):
        if num_cb == 0:
            return gate

        q_2 = gate.mat.shape[0]
        q_1 = q_2 // 2
        q_3 = q_2 * 3 // 2
        q_4 = q_2 * 2

        quads = [gate.mat[:q_1, :q_1], gate.mat[q_1:, :q_1], gate.mat[:q_1, q_1:], gate.mat[q_1:, q_1:]]

        gate.mat = np.zeros(shape=(q_4,q_4))
        gate.mat[0:q_1, 0:q_1] += np.eye(q_1)
        gate.mat[q_2:q_3, q_2:q_3] += np.eye(q_1)

        gate.mat[q_1:q_2, q_1:q_2] += quads[0]
        gate.mat[q_3:q_4, q_1:q_2] += quads[1]
        gate.mat[q_1:q_2, q_3:q_4] += quads[2]
        gate.mat[q_3:q_4, q_3:q_4] += quads[3]

        gate.num_cb += 1
        return Gates.C_mat(gate, num_cb - 1)

    @staticmethod
    def check_size(mat, i_qbits):
        len_i_qbits = len(i_qbits)
        len_matrix = np.log2(mat.shape[0])

        if len_i_qbits != len_matrix:
            print("Error: Dimension mismatch: " + str(len_i_qbits) + " " + str(len_matrix))
            return False
        return True

    @staticmethod
    def check_set(i_qbits):
        if len(i_qbits) is not len(list(dict.fromkeys(i_qbits))):
            print("Error: i_qbits is not a set: " + i_qbits.__repr__())
            return False

        return True

    @staticmethod
    def check_valid_i_qbits(i_qbits, n_qbits):
        for i_qbit in i_qbits:
            if i_qbit >= n_qbits:
                print("Error: Invalid i_qbits values: " + i_qbits.__repr__())
                return False

        return True

    @staticmethod
    def print_gate(gate_inst, i_qbit):
        num_chars = len(gate_inst.gate.gate_label)

        if i_qbit not in gate_inst.i_qbits:
            [print(" ", end="") for i in range(num_chars)]
            if gate_inst.gate.num_cb > 0:
                print(" ", end="")
            return

        j_qbit = gate_inst.i_qbits.index(i_qbit)
        if j_qbit < gate_inst.gate.num_cb:
            print("C", end="")
            [print(" ", end="") for i in range(num_chars)]
        else:
            if gate_inst.gate.num_cb > 0:
                print(" ", end="")
            print(gate_inst.gate.gate_label, end="")