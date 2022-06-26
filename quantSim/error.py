import numpy as np


class Error:

    @staticmethod
    def check_size(mat, i_qbits):
        """
        Ensures the size of a matrix corresponds to the correct number of qubits
        """
        len_i_qbits = len(i_qbits)
        len_matrix = np.log2(mat.shape[0])

        if len_i_qbits != len_matrix:
            print(f"Error: Dimension mismatch: {str(len_i_qbits)} {str(len_matrix)}")
            return False
        return True

    @staticmethod
    def check_set(i_qbits):
        """
        Ensures i_qbits is a set
        """
        if len(i_qbits) is not len(list(dict.fromkeys(i_qbits))):
            print(f"Error: i_qbits is not a set: {i_qbits.__repr__()}")
            return False

        return True

    @staticmethod
    def check_valid_i_qbits(i_qbits, n_qbits):
        """
        Ensures indices in i_qbits are valid (lower than n_qbits)
        """
        for i_qbit in i_qbits:
            if i_qbit >= n_qbits:
                print(f"Error: Invalid i_qbits values: {i_qbits.__repr__()}")
                return False

        return True
