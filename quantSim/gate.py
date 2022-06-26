import numpy as np
from quantSim import GateContainer


class Gate:
    # common constants
    SQRT_2 = 1 / np.sqrt(2)
    SQRT_8 = 1 / np.sqrt(8)
    SQRT_32 = 1 / np.sqrt(32)

    @staticmethod
    def I():
        """
        :return: an Identity gate
        """
        mat = np.eye(2)
        return GateContainer("I", mat, 0)

    @staticmethod
    def H():
        """
        :return: a Hadamard gate
        """
        mat = np.array(
            [[Gate.SQRT_2, Gate.SQRT_2],
             [Gate.SQRT_2, -Gate.SQRT_2]])

        return GateContainer("H", mat, 0)

    @staticmethod
    def X():
        """
        :return: a not gate
        """
        mat = np.array(
            [[0, 1],
             [1, 0]])

        return GateContainer("X", mat, 0)

    @staticmethod
    def Y():
        """
        :return: a Y gate
        """
        mat = np.array(
            [[0, -1j],
             [1j, 0]])

        return GateContainer("Y", mat, 0)

    @staticmethod
    def Z():
        """
        :return: a Z gate
        """
        mat = np.array(
            [[1, 0],
             [0, -1]])

        return GateContainer("Z", mat, 0)

    @staticmethod
    def rX():
        """
        :return: a square root of not gate
        """
        mat = np.array(
            [[0.5 + 0.5j, 0.5 - 0.5j],
             [0.5 - 0.5j, 0.5 + 0.5j]])

        return GateContainer("rX", mat, 0)

    @staticmethod
    def P(theta, gate_label="P"):
        """
        :param theta: phase shift angle in radians
        :return: a phase shift gate
        """
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

    @staticmethod
    def S(dagger=False):
        """
        :return: an S gate
        """
        if dagger:
            return Gate.Sdg()

        mat = np.array(
            [[1, 0],
             [0, 1j]])

        return GateContainer("S", mat, 0)

    @staticmethod
    def Sdg():
        """
        :return: an S dagger gate
        """
        mat = np.array(
            [[1, 0],
             [0, -1j]])

        return GateContainer("Sdg", mat, 0)

    @staticmethod
    def T(dagger=False):
        """
        :return: a T gate
        """
        if dagger:
            return Gate.Tdg()

        mat = np.array(
            [[1, 0],
             [0, Gate.SQRT_2 + 1j * Gate.SQRT_2]])

        return GateContainer("T", mat, 0)

    @staticmethod
    def Tdg():
        """
        :return: a T dagger gate
        """
        mat = np.array(
            [[1, 0],
             [0, Gate.SQRT_2 - 1j * Gate.SQRT_2]])

        return GateContainer("Tdg", mat, 0)

    @staticmethod
    def Swap():
        """
        :return: a swap gate
        """
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [0, 0, 0, 1]])

        return GateContainer("SW", mat, 0)

    @staticmethod
    def Swap_sqrt():
        """
        :return: a square root of swap gate
        """
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0.5 + 0.5j, 0.5 - 0.5j, 0],
             [0, 0.5 - 0.5j, 0.5 + 0.5j, 0],
             [0, 0, 0, 1]])

        return GateContainer("rSW", mat, 0)

    @staticmethod
    def Swap_i():
        """
        :return: an imaginary swap gate
        """
        mat = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1j, 0],
             [0, 1j, 0, 0],
             [0, 0, 0, 1]])

        return GateContainer("iSW", mat, 0)

    @staticmethod
    def Swap_sqrt_i():
        """
        :return: an imaginary square root of swap gate
        """
        mat = np.array(
            [[1, 0, 0, 0],
             [0, Gate.SQRT_2, 1j * Gate.SQRT_2, 0],
             [0, 1j * Gate.SQRT_2, Gate.SQRT_2, 0],
             [0, 0, 0, 1]])

        return GateContainer("irSW", mat, 0)

    """
        Gates created with other matrices
    """

    @staticmethod
    def CX():
        """
        :return: a controlled not gate
        """
        gate = Gate.X()
        return Gate.MC(gate, num_cb=1)

    @staticmethod
    def Toffoli():
        """
        :return: a Toffoli gate
        """
        gate = Gate.X()
        return Gate.MC(gate, num_cb=2)

    @staticmethod
    def CSwap():
        """
        :return: a controlled swap gate
        """
        gate = Gate.Swap()
        return Gate.MC(gate, num_cb=1)

    @staticmethod
    def CP(theta):
        """
        :param theta: phase shift angle in radians
        :return: a controlled Phase gate
        """
        gate = Gate.P(theta)
        return Gate.MC(gate, num_cb=1)

    """
        Special gate operations
    """

    @staticmethod
    def dg(gate):
        """
        :param gate: modification of this
        :return: conjugate transpose of gate
        """
        mat = np.matrix.h(gate.mat)

        return GateContainer(gate.gate_label, mat, gate.num_cb)


    @staticmethod
    def MT(gate, num_tb):
        """
        :param gate: target gate
        :param num_tb: number of target gates
        :return: a multi-target gate version of a given single bit gate
        """
        n_qbits_mat = np.log2(gate.mat.shape[0])
        if n_qbits_mat != 1:
            print(f"Error: Gate is not single bit: {n_qbits_mat}")
            return

        mat = [1]
        for i in range(num_tb):
            mat = np.kron(mat, gate.mat)

        return GateContainer(gate.gate_label, mat, 0)

    @staticmethod
    def MC(gate, num_cb):
        """
        :param gate: target gate
        :param num_cb: number of control bits
        :return: a controlled version of gate recursively
        """
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