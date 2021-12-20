import numpy as np
#import sim

"""
    gate_label: characters to be printed when circuit diagram is printed (not including the 'C's for control bits)
    mat: gate matrix
        -LSB
        -includes control bits
    num_cb: number of control bits for this gate, control bits are the MSBs, only used for printing
    num_args: number of arguments the gate is expecting
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

    def __init__(self, gate, i_qbits, n_qbits):
        self.gate = gate
        self.i_qbits = i_qbits
        self.n_qbits = n_qbits

        if Gates.check_size(gate.mat, i_qbits):
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

    # Gate labels
    # 1 qubit
    I = "I"
    X = "X"
    Y = "Y"
    Z = "Z"
    H = "H"
    S = "S"
    T = "T"
    # 2 qubit
    CX = "CX"
    CY = "CY"
    CZ = "CZ"
    SWAP = "SW"
    sqSWAP = "sqSW"
    iSWAP = "iSW"
    sqiSWAP = "sqiSW"
    # 3 qubit
    TOFFOLI = "Tof"
    CSWAP = "CSW"

    # Gate matrices
    GATE_I_M = np.array(
        [[1, 0],
         [0, 1]])

    GATE_X_M = np.array(
        [[0, 1],
         [1, 0]])

    GATE_Y_M = np.array(
        [[0, -1j],
         [1j, 0]])

    GATE_Z_M = np.array(
        [[1, 0],
         [0, -1]])

    GATE_H_M = np.array(
        [[SQRT_H, SQRT_H],
         [SQRT_H, -SQRT_H]])

    GATE_S_M = np.array(
        [[1, 0],
         [0, 1j]])

    GATE_T_M = np.array(
        [[1, 0],
         [0, pow(np.e, np.pi * 1j / 4)]])

    GATE_CX_M = np.array(
        [[1, 0, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 1, 0],
         [0, 1, 0, 0]])

    GATE_CY_M = np.array(
        [[1, 0, 0, 0],
         [0, 0, 0, -1j],
         [0, 0, 1, 0],
         [0, 1j, 0, 0]])

    GATE_CZ_M = np.array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, -1]])

    GATE_SWAP_M = np.array(
        [[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1]],dtype=complex)

    GATE_sqSWAP_M = np.array(
        [[1, 0, 0, 0],
         [0, 0.5 + 0.5j, 0.5 - 0.5j, 0],
         [0, 0.5 - 0.5j, 0.5 + 0.5j, 0],
         [0, 0, 0, 1]])

    GATE_iSWAP_M = np.array(
        [[1, 0, 0, 0],
         [0, 0, 1j, 0],
         [0, 1j, 0, 0],
         [0, 0, 0, 1]])

    GATE_sqiSWAP_M = np.array(
        [[1, 0, 0, 0],
         [0, SQRT_H, SQRT_H * 1j, 0],
         [0, SQRT_H * 1j, SQRT_H, 0],
         [0, 0, 0, 1]])

    GATE_TOFFOLI_M = np.array(
        [[1, 0, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 1, 0, 0, 0, 0]])

    GATE_CSWAP_M = np.array(
        [[1, 0, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 1]])

    gate_dict = {}
    gate_dict[I] = GATE_I_M
    gate_dict[X] = GATE_X_M
    gate_dict[Y] = GATE_Y_M
    gate_dict[Z] = GATE_Z_M
    gate_dict[H] = GATE_H_M
    gate_dict[S] = GATE_S_M
    gate_dict[T] = GATE_T_M
    gate_dict[CX] = GATE_CX_M
    gate_dict[CY] = GATE_CY_M
    gate_dict[CZ] = GATE_CZ_M
    gate_dict[SWAP] = GATE_SWAP_M
    gate_dict[sqSWAP] = GATE_sqSWAP_M
    gate_dict[iSWAP] = GATE_iSWAP_M
    gate_dict[sqiSWAP] = GATE_sqiSWAP_M
    gate_dict[TOFFOLI] = GATE_TOFFOLI_M
    gate_dict[CSWAP] = GATE_CSWAP_M

    # number of character to print
    gate_print = {}
    gate_print[I] = 1
    gate_print[X] = 1
    gate_print[Y] = 1
    gate_print[Z] = 1
    gate_print[H] = 1
    gate_print[S] = 1
    gate_print[T] = 1
    gate_print[CX] = 2
    gate_print[CY] = 2
    gate_print[CZ] = 2
    gate_print[SWAP] = 2
    gate_print[sqSWAP] = 4
    gate_print[iSWAP] = 3
    gate_print[sqiSWAP] = 5
    gate_print[TOFFOLI] = 2
    gate_print[CSWAP] = 3

    # which gates have control/controlled qubits
    # all qubits are control except last
    gate_print_control = {}
    gate_print_control[CX] = "X"
    gate_print_control[CY] = "Y"
    gate_print_control[CZ] = "Z"
    gate_print_control[TOFFOLI] = "X"
    gate_print_control[CSWAP] = "SW"

    @staticmethod
    def GATE_I(i_qbits, n_qbits):
        gate = Gates.I_mat()

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_H(i_qbits, n_qbits):
        gate = Gates.H_mat()

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_X(i_qbits, n_qbits):
        gate = Gates.X_mat()

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_Y(i_qbits, n_qbits):
        gate = Gates.Y_mat()

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_Z(i_qbits, n_qbits):
        gate = Gates.Z_mat()

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_CX(i_qbits, n_qbits):
        num_cb = 1
        gate = Gates.X_mat()
        gate = Gates.C_mat(gate, num_cb)

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_TOFFOLI(i_qbits, n_qbits):
        num_cb = 2
        gate = Gates.X_mat()
        gate = Gates.C_mat(gate, num_cb)

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_P(i_qbits, n_qbits, theta):
        gate = Gates.P_mat(theta)

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)

    @staticmethod
    def GATE_CP(i_qbits, n_qbits, theta):
        num_cb = 1
        gate = Gates.P_mat(theta)
        gate = Gates.C_mat(gate, num_cb)

        if Gates.check_size(gate.mat, i_qbits):
            return Gate_Inst(gate, i_qbits, n_qbits)


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
    def T_mat():
        return Gates.P_mat(np.pi / 4, "T")

    @staticmethod
    def P_mat(theta, gate_label="P"):
        mat = np.array(
            [[1, 0],
             [0, np.exp(theta * 1j)]])

        return Gate(gate_label, mat, 0)

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

    # creates a controlled version of base_mat with num_cb number of control bits
    @staticmethod
    def C_mat(gate, num_cb):
        base_n_qbits = np.log2(gate.mat.shape[0])

        return gate

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
    def print_gate_2(gate_inst, i_qbit):
        num_chars = len(gate_inst.gate.gate_label)

        if i_qbit not in i_qbit:
            [print(" ", end="") for i in range(num_chars)]
            if gate_inst.gate.num_cb > 0:
                print(" ", end="")
            return

        j_qbit = gate_inst.i_qbits.index(i_qbit)
        if j_qbit < gate_inst.gate.num_cb:
            print("C", end="")
            [print(" ", end="") for i in range(num_chars)]
        else:
            print(" " + gate_inst.gate.gate_label, end="")


    @staticmethod
    def print_gate(gate_label, i_qbits, i_qbit):
        num_chars = Gates.gate_print[gate_label]

        if i_qbit not in i_qbits:
            [print(" ", end="") for i in range(num_chars)]
            return

        if num_chars == 1:
            print(f'{gate_label:1}', end="")
        elif num_chars == 2:
            if gate_label in Gates.gate_print_control:
                if i_qbit == i_qbits[len(i_qbits) - 1]:
                    print(f' {Gates.gate_print_control[gate_label]:1}', end="")
                else:
                    print("C ", end="")
            else:
                print(gate_label, end="")
        else:
            print(gate_label, end="")