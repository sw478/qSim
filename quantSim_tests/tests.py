import numpy as np
import sys
import pytest
import pytest_runner
from quantSim import *

"""
    Some gates are hardcoded in Gate but still have tests since project architecture is subject to change
"""

margin_of_error = 0.0001

"""
    Checking gate matrices of some single-bit gates:
    - I, H, X, Y, Z, S, Sdg, T, Tdg
    - P: has an argument theta:
        - checking equivalence to Z, S, Sdg, T, Tdg
"""
def test_gates_basic():

    expected_I = np.array(
        [[1, 0],
         [0, 1]])

    expected_H = np.array(
        [[Gate.SQRT_2, Gate.SQRT_2],
         [Gate.SQRT_2, -Gate.SQRT_2]])

    expected_X = np.array(
        [[0, 1],
         [1, 0]])

    expected_Y = np.array(
        [[0, -1j],
         [1j, 0]])

    expected_Z = np.array(
        [[1, 0],
         [0, -1]])

    expected_S = np.array(
        [[1, 0],
         [0, 1j]])

    expected_Sdg = np.array(
        [[1, 0],
         [0, -1j]])

    expected_T = np.array(
        [[1, 0],
         [0, Gate.SQRT_2 + 1j * Gate.SQRT_2]])

    expected_Tdg = np.array(
        [[1, 0],
         [0, Gate.SQRT_2 - 1j * Gate.SQRT_2]])

    assert_matrices(Gate.I().mat, expected_I)
    assert_matrices(Gate.H().mat, expected_H)
    assert_matrices(Gate.X().mat, expected_X)
    assert_matrices(Gate.Y().mat, expected_Y)
    assert_matrices(Gate.Z().mat, expected_Z)

    assert_matrices(Gate.S().mat, expected_S)
    assert_matrices(Gate.T().mat, expected_T)
    assert_matrices(Gate.Sdg().mat, expected_Sdg)
    assert_matrices(Gate.Tdg().mat, expected_Tdg)
    assert_matrices(Gate.S(True).mat, expected_Sdg)
    assert_matrices(Gate.T(True).mat, expected_Tdg)

    theta_Z = np.pi
    theta_S = np.pi / 2
    theta_T = np.pi / 4

    assert_matrices(Gate.P(theta_Z).mat, expected_Z)
    assert_matrices(Gate.P(theta_S).mat, expected_S)
    assert_matrices(Gate.P(-theta_S).mat, expected_Sdg)
    assert_matrices(Gate.P(theta_T).mat, expected_T)
    assert_matrices(Gate.P(-theta_T).mat, expected_Tdg)


"""
    Checking gate matrices of rotational gates:
    - Rx, Ry, Rz
    - Gates have an argument theta:
        - 0, pi/2, pi, 2pi, 4pi
"""


def test_rotational_gates():
    theta_0 = 0
    theta_1 = np.pi / 2
    theta_2 = np.pi
    theta_4 = np.pi * 2
    theta_8 = np.pi * 4

    expected_Rx_0 = np.array(
        [[1, 0],
         [0, 1]])

    expected_Rx_1 = np.array(
        [[Gate.SQRT_2, -1j * Gate.SQRT_2],
         [-1j * Gate.SQRT_2, Gate.SQRT_2]])

    expected_Rx_2 = np.array(
        [[0, -1j],
         [-1j, 0]])

    expected_Rx_4 = np.array(
        [[-1, 0],
         [0, -1]])

    assert_matrices(Gate.Rx(theta_0).mat, expected_Rx_0)
    assert_matrices(Gate.Rx(theta_1).mat, expected_Rx_1)
    assert_matrices(Gate.Rx(theta_2).mat, expected_Rx_2)
    assert_matrices(Gate.Rx(theta_4).mat, expected_Rx_4)
    assert_matrices(Gate.Rx(theta_8).mat, expected_Rx_0)

    expected_Ry_0 = np.array(
        [[1, 0],
         [0, 1]])

    expected_Ry_1 = np.array(
        [[Gate.SQRT_2, -Gate.SQRT_2],
         [Gate.SQRT_2, Gate.SQRT_2]])

    expected_Ry_2 = np.array(
        [[0, -1],
         [1, 0]])

    expected_Ry_4 = np.array(
        [[-1, 0],
         [0, -1]])

    assert_matrices(Gate.Ry(theta_0).mat, expected_Ry_0)
    assert_matrices(Gate.Ry(theta_1).mat, expected_Ry_1)
    assert_matrices(Gate.Ry(theta_2).mat, expected_Ry_2)
    assert_matrices(Gate.Ry(theta_4).mat, expected_Ry_4)
    assert_matrices(Gate.Ry(theta_8).mat, expected_Ry_0)

    expected_Rz_0 = np.array(
        [[1, 0],
         [0, 1]])

    expected_Rz_1 = np.array(
        [[Gate.SQRT_2 - 1j * Gate.SQRT_2, 0],
         [0, Gate.SQRT_2 + 1j * Gate.SQRT_2]])

    expected_Rz_2 = np.array(
        [[-1j, 0],
         [0, 1j]])

    expected_Rz_4 = np.array(
        [[-1, 0],
         [0, -1]])

    assert_matrices(Gate.Rz(theta_0).mat, expected_Rz_0)
    assert_matrices(Gate.Rz(theta_1).mat, expected_Rz_1)
    assert_matrices(Gate.Rz(theta_2).mat, expected_Rz_2)
    assert_matrices(Gate.Rz(theta_4).mat, expected_Rz_4)
    assert_matrices(Gate.Rz(theta_8).mat, expected_Rz_0)


"""
    Checking gate matrices of ising coupling gates:
    - Rxx, Ryy, Rzz
    - Gates have an argument theta:
        - 0, pi/2, pi, 2pi, 4pi
"""


def test_ising_coupling_gates():
    # increments of PI/2
    theta_0 = 0
    theta_1 = np.pi / 2
    theta_2 = np.pi
    theta_4 = np.pi * 2
    theta_8 = np.pi * 4

    diag0 = 1
    diag1 = 0
    expected_Rxx_0 = np.array(
        [[diag0, 0, 0, diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [diag1, 0, 0, diag0]])

    diag0 = Gate.SQRT_2
    diag1 = -1j * Gate.SQRT_2
    expected_Rxx_1 = np.array(
        [[diag0, 0, 0, diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [diag1, 0, 0, diag0]])

    diag0 = 0
    diag1 = -1j
    expected_Rxx_2 = np.array(
        [[diag0, 0, 0, diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [diag1, 0, 0, diag0]])

    diag0 = -1
    diag1 = 0
    expected_Rxx_4 = np.array(
        [[diag0, 0, 0, diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [diag1, 0, 0, diag0]])

    assert_matrices(Gate.Rxx(theta_0).mat, expected_Rxx_0)
    assert_matrices(Gate.Rxx(theta_1).mat, expected_Rxx_1)
    assert_matrices(Gate.Rxx(theta_2).mat, expected_Rxx_2)
    assert_matrices(Gate.Rxx(theta_4).mat, expected_Rxx_4)
    assert_matrices(Gate.Rxx(theta_8).mat, expected_Rxx_0)

    diag0 = 1
    diag1 = 0
    expected_Ryy_0 = np.array(
        [[diag0, 0, 0, -diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [-diag1, 0, 0, diag0]])

    diag0 = Gate.SQRT_2
    diag1 = -1j * Gate.SQRT_2
    expected_Ryy_1 = np.array(
        [[diag0, 0, 0, -diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [-diag1, 0, 0, diag0]])

    diag0 = 0
    diag1 = -1j
    expected_Ryy_2 = np.array(
        [[diag0, 0, 0, -diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [-diag1, 0, 0, diag0]])

    diag0 = -1
    diag1 = 0
    expected_Ryy_4 = np.array(
        [[diag0, 0, 0, -diag1],
         [0, diag0, diag1, 0],
         [0, diag1, diag0, 0],
         [-diag1, 0, 0, diag0]])

    assert_matrices(Gate.Ryy(theta_0).mat, expected_Ryy_0)
    assert_matrices(Gate.Ryy(theta_1).mat, expected_Ryy_1)
    assert_matrices(Gate.Ryy(theta_2).mat, expected_Ryy_2)
    assert_matrices(Gate.Ryy(theta_4).mat, expected_Ryy_4)
    assert_matrices(Gate.Ryy(theta_8).mat, expected_Ryy_0)

    diag0 = 1
    diag1 = 1
    expected_Rzz_0 = np.array(
        [[diag0, 0, 0, 0],
         [0, diag1, 0, 0],
         [0, 0, diag1, 0],
         [0, 0, 0, diag0]])

    diag0 = Gate.SQRT_2 - 1j * Gate.SQRT_2
    diag1 = Gate.SQRT_2 + 1j * Gate.SQRT_2
    expected_Rzz_1 = np.array(
        [[diag0, 0, 0, 0],
         [0, diag1, 0, 0],
         [0, 0, diag1, 0],
         [0, 0, 0, diag0]])

    diag0 = -1j
    diag1 = 1j
    expected_Rzz_2 = np.array(
        [[diag0, 0, 0, 0],
         [0, diag1, 0, 0],
         [0, 0, diag1, 0],
         [0, 0, 0, diag0]])

    diag0 = -1
    diag1 = -1
    expected_Rzz_4 = np.array(
        [[diag0, 0, 0, 0],
         [0, diag1, 0, 0],
         [0, 0, diag1, 0],
         [0, 0, 0, diag0]])

    assert_matrices(Gate.Rzz(theta_0).mat, expected_Rzz_0)
    assert_matrices(Gate.Rzz(theta_1).mat, expected_Rzz_1)
    assert_matrices(Gate.Rzz(theta_2).mat, expected_Rzz_2)
    assert_matrices(Gate.Rzz(theta_4).mat, expected_Rzz_4)
    assert_matrices(Gate.Rzz(theta_8).mat, expected_Rzz_0)


"""
    Checking gate matrices of swap gates:
    - Swap, sqrt Swap, imag Swap, sqrt imag Swap
"""


def test_swap_and_variations():
    expected_swap = np.array(
        [[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1]])

    expected_swap_sqrt = np.array(
        [[1, 0, 0, 0],
         [0, 0.5 + 0.5j, 0.5 - 0.5j, 0],
         [0, 0.5 - 0.5j, 0.5 + 0.5j, 0],
         [0, 0, 0, 1]])

    expected_swap_i = np.array(
        [[1, 0, 0, 0],
         [0, 0, 1j, 0],
         [0, 1j, 0, 0],
         [0, 0, 0, 1]])

    expected_swap_sqrt_i = np.array(
        [[1, 0, 0, 0],
         [0, Gate.SQRT_2, 1j * Gate.SQRT_2, 0],
         [0, 1j * Gate.SQRT_2, Gate.SQRT_2, 0],
         [0, 0, 0, 1]])

    assert_matrices(Gate.Swap().mat, expected_swap)
    assert_matrices(Gate.Swap_sqrt().mat, expected_swap_sqrt)
    assert_matrices(Gate.Swap_i().mat, expected_swap_i)
    assert_matrices(Gate.Swap_sqrt_i().mat, expected_swap_sqrt_i)


"""
    Checking gate matrices of control gates:
    - CX, Toffoli, CSwap, CPhase
"""


def test_control_gates():
    expected_cx = np.array(
        [[1, 0, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 1, 0],
         [0, 1, 0, 0]])

    expected_toffoli = np.array(
        [[1, 0, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 1, 0, 0, 0, 0]])

    expected_cswap = np.array(
        [[1, 0, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 1]])

    assert_matrices(Gate.CX().mat, expected_cx)
    assert_matrices(Gate.Toffoli().mat, expected_toffoli)
    assert_matrices(Gate.CSwap().mat, expected_cswap)

    expected_cp_z = np.array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, -1]])

    expected_cp_s = np.array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, 1j]])

    expected_cp_t = np.array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, Gate.SQRT_2 + 1j * Gate.SQRT_2]])

    theta_Z = np.pi
    theta_S = np.pi / 2
    theta_T = np.pi / 4

    assert_matrices(Gate.CP(theta_Z).mat, expected_cp_z)
    assert_matrices(Gate.CP(theta_S).mat, expected_cp_s)
    assert_matrices(Gate.CP(theta_T).mat, expected_cp_t)


"""
    Checking statevector of Sim with multi-target gates
"""

def test_multi_target_gates():
    n_qbits = 4
    sim = Sim(n_qbits, "Test: Multi-Target Gates", True)
    print("")

    deg_045 = 0.25 * np.exp(1j * 1 * np.pi / 4)
    deg_135 = 0.25 * np.exp(1j * 3 * np.pi / 4)
    deg_225 = 0.25 * np.exp(1j * 5 * np.pi / 4)
    deg_315 = 0.25 * np.exp(1j * 7 * np.pi / 4)
    expected_prob_dist = [0.0625 for i in range(1 << n_qbits)]

    sim.add_gate(Gate.MT(Gate.H(), 4), [0, 1, 2, 3])
    assert_arrays(sim.print_prob_dist(), expected_prob_dist)
    assert_arrays(sim.print_statevector(), [0.25 for i in range(1 << n_qbits)])

    sim.add_gate(Gate.MT(Gate.T(), 4), [0, 1, 2, 3])
    assert_arrays(sim.print_prob_dist(), expected_prob_dist)
    assert_arrays(sim.print_statevector(), [0.25, deg_045, deg_045, 0.25j, deg_045, 0.25j, 0.25j, deg_135, deg_045, 0.25j, 0.25j, deg_135, 0.25j, deg_135, deg_135, -0.25])

    sim.add_gate(Gate.MT(Gate.X(), 3), [0, 1, 3])
    assert_arrays(sim.print_prob_dist(), expected_prob_dist)
    assert_arrays(sim.print_statevector(), [deg_135, 0.25j, 0.25j, deg_045, -0.25, deg_135, deg_135, 0.25j, 0.25j, deg_045, deg_045, 0.25, deg_135, 0.25j, 0.25j, deg_045])

    sim.add_gate(Gate.MT(Gate.S(), 3), [0, 2, 3])
    assert_arrays(sim.print_prob_dist(), expected_prob_dist)
    assert_arrays(sim.print_statevector(), [deg_135, -0.25, 0.25j, deg_135, -0.25j, deg_315, deg_225, -0.25j, -0.25, deg_225, deg_135, -0.25, deg_315, 0.25, -0.25j, deg_315])

    sim.add_gate(Gate.MT(Gate.Y(), 2), [1, 2])
    assert_arrays(sim.print_prob_dist(), expected_prob_dist)
    assert_arrays(sim.print_statevector(), [deg_045, 0.25j, -0.25j, deg_315, 0.25j, deg_135, deg_315, 0.25, 0.25j, deg_135, deg_315, 0.25, deg_135, -0.25, 0.25, deg_045])

    sim.print_sim()


"""
    Checking full matrices of SubGate
"""


def test_sub_gates():
    print("")
    n_qbits = 3

    sub_gate = SubGate(Gate.X(), [2], n_qbits)

    base_mat = sub_gate.gate.mat
    expected = np.array([[0, 1], [1, 0]])
    assert_matrices(base_mat, expected)

    full_mat = sub_gate.full_mat
    expected = np.array(
        [[0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 1],
         [1, 0, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0]])
    assert_matrices(full_mat, expected)


"""
    Checking statevector of Sim with gates
"""


def test_statevector():
    n_qbits = 3
    sim = Sim(n_qbits, "Test: Statevector", True)
    print("")

    gate = Gate.MT(Gate.H(), 2)
    sim.add_gate(gate, [0, 1])
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0])

    sim.add_gate(Gate.Toffoli(), [0, 1, 2])
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0, 0, 0, 0, 0.5])

    sim.add_gate(Gate.Swap(), [1, 2])
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0, 0, 0.5, 0, 0, 0.5])

    sim.add_gate(Gate.H(), [2])
    assert_arrays(sim.print_statevector(), [Gate.SQRT_2, Gate.SQRT_8, 0, Gate.SQRT_8, 0, Gate.SQRT_8, 0, -Gate.SQRT_8])

    sim.add_gate(Gate.Z(), [0])
    assert_arrays(sim.print_statevector(),
                  [Gate.SQRT_2, -Gate.SQRT_8, 0, -Gate.SQRT_8, 0, -Gate.SQRT_8, 0, Gate.SQRT_8])

    sim.add_gate(Gate.S(), [1])
    assert_arrays(sim.print_statevector(),
                  [Gate.SQRT_2, -Gate.SQRT_8, 0, -1j * Gate.SQRT_8, 0, -Gate.SQRT_8, 0, 1j * Gate.SQRT_8])

    sim.add_gate(Gate.T(), [2])
    assert_arrays(sim.print_statevector(),
                  [Gate.SQRT_2, -Gate.SQRT_8, 0, -Gate.SQRT_8 * 1j, 0, -0.25 - 0.25j, 0, -0.25 + 0.25j])

    sim.print_sim()


"""
    Checking statevector of custom Gates made from Sims
"""


def test_custom_gate():
    n_qbits_1 = 3
    sim1 = Sim(n_qbits_1, "Custom Gate", True)
    print("")

    sim1.add_gate(Gate.H(), [0])
    sim1.add_gate(Gate.H(), [1])
    sim1.add_gate(Gate.Toffoli(), [0, 1, 2])
    sim1.add_gate(Gate.Swap(), [1, 2])

    sim1.print_sim()
    gate_sim1 = sim1.to_gate("G1")

    assert gate_sim1.expected_qbits == 3

    n_qbits_2 = 7
    sim2 = Sim(n_qbits_2, "Custom Gate")
    print("")

    sim2.add_gate(gate_sim1, [0, 1, 2])
    sim2.add_gate(Gate.Toffoli(), [1, 2, 4])
    sim2.add_gate(Gate.Swap(), [1, 4])
    sim2.add_gate(Gate.Swap(), [2, 6])
    sim2.add_gate(Gate.Toffoli(), [2, 5, 3])

    sim2.print_sim()
    # sim2.print_statevector()


def assert_arrays(arr1, arr2):
    assert arr1 == pytest.approx(arr2, margin_of_error)


def assert_matrices(a, b):
    np.testing.assert_array_almost_equal(a, b, 4)


def main():
    pass
    # running with pytest


if __name__ == "__main__":
    np.set_printoptions(threshold=sys.maxsize)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
    main()
