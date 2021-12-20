import numpy as np
from inspect import signature

"""
    gate_label: characters to be printed when circuit diagram is printed (not including the 'C's for control bits)
    mat: function pointer for matrix
        -arguments for function are arguments for matrix
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
        self.num_args = len(signature(mat).parameters)

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
    def X():
        return np.array(
            [[0, 1],
             [1, 0]])

    @staticmethod
    def PHASE(theta):
        return np.array(
        [[1, 0],
         [0, pow(np.e, theta * 1j)]])

    @staticmethod
    def Rx(theta):
        theta_2 = theta * 0.5
        return np.array(
            [[np.cos(theta_2), np.sin(theta_2) * -1j],
             [np.sin(theta_2) * -1j, np.cos(theta_2)]]
        )

    @staticmethod
    def Ry(theta):
        theta_2 = theta * 0.5
        return np.array(
            [[np.cos(theta_2), np.sin(theta_2) * -1],
             [np.sin(theta_2), np.cos(theta_2)]]
        )

    @staticmethod
    def Rz(theta):
        theta_2 = theta * 0.5
        return np.array(
            [[np.exp(-1j * theta_2), 0],
             [0, np.exp(1j * theta_2)]]
        )

    @staticmethod
    def Rxx(theta):
        theta_2 = theta * 0.5
        cos_theta_2 = np.cos(theta_2)
        isin_theta_2_neg = -1j * np.sin(theta_2)

        return np.array(
            [[cos_theta_2, 0, 0, isin_theta_2_neg],
             [0, cos_theta_2, isin_theta_2_neg, 0],
             [0, isin_theta_2_neg, cos_theta_2, 0],
             [isin_theta_2_neg, 0, 0, cos_theta_2]]
        )

    @staticmethod
    def Ryy(theta):
        theta_2 = theta * 0.5
        cos_theta_2 = np.cos(theta_2)
        isin_theta_2_neg = -1j * np.sin(theta_2)

        return np.array(
            [[cos_theta_2, 0, 0, -isin_theta_2_neg],
             [0, cos_theta_2, isin_theta_2_neg, 0],
             [0, isin_theta_2_neg, cos_theta_2, 0],
             [-isin_theta_2_neg, 0, 0, cos_theta_2]]
        )

    @staticmethod
    def Rzz(theta):
        theta_2 = theta * 0.5
        e_theta_2 = np.exp(1j * theta_2)
        e_theta_2_neg = np.exp(-1j * theta_2)

        return np.array(
            [[e_theta_2_neg, 0, 0, 0],
             [0, e_theta_2, 0, 0],
             [0, 0, e_theta_2, 0],
             [0, 0, 0, e_theta_2_neg]]
        )

    # creates a controlled version of mat with num_control_bits number of control bits
    @staticmethod
    def C(mat, num_control_bits):
        pass

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