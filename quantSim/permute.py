import numpy as np


class Permute:
    __swap_mat = np.array(
        [[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1]])

    __max_n_qbits = 0
    __preproc = None

    @staticmethod
    def get_swap_mat(n_qbits, i_qbit, j_qbit):
        """
        :param n_qbits: extended number of qubits
        :return: a swap matrix to swap qubits (i_qbit and j_qbit)
        """
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
        """
        Preprocesses all swap matrices for indices between 0 and n_qbits - 1

        :param n_qbits: max number of qubits to preprocess for now
        """
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
                    mat = np.kron(Permute.__swap_mat, mat)
                else:
                    mat = np.kron(np.eye(2), mat)

            Permute.__preproc[i][i + 1] = mat
            Permute.__preproc[i + 1][i] = mat

        for j in range(2, n_qbits):
            for i in range(j - 2, -1, -1):
                a = Permute.__preproc[i][i + 1]
                b = Permute.__preproc[i + 1][j]
                Permute.__preproc[i][j] = np.linalg.multi_dot([a, b, a])

