import numpy as np
import sys
import pytest
from gates import Gates
from sim import Sim

margin = 0.000000001

def test_gates():
    print("")
    n_qbits = 3
    i_qbits = [2]
    gate_inst = Gates.X(i_qbits, n_qbits)

    base_mat = gate_inst.gate.mat
    expected = np.array([[0, 1], [1, 0]])
    assert_matrices(base_mat, expected)

    full_mat = gate_inst.full_mat
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

def test_gates_with_args():
    print("")
    n_qbits = 3
    i_qbits = [2]
    theta = np.pi / 2
    gate_inst = Gates.P(i_qbits, n_qbits, theta)

    base_mat = gate_inst.gate.mat
    expected = np.array([[1, 0], [0, 1j]])
    assert_matrices(base_mat, expected)

def test_controlled_gates():
    n_qbits = 3
    sim = Sim(n_qbits, "Test: Controlled Gates", True)
    print("")

    sim.add_gate(Gates.H([0], n_qbits))
    sim.add_gate(Gates.H([1], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0])

    sim.add_gate(Gates.Toffoli([0, 1, 2], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0, 0, 0, 0, 0.5])


def test_statevector():
    n_qbits = 3
    sim = Sim(n_qbits, "Test: Statevector", True)
    print("")

    sim.add_gate(Gates.H([0], n_qbits))
    sim.add_gate(Gates.H([1], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0])

    sim.add_gate(Gates.Toffoli([0, 1, 2], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0, 0, 0, 0, 0.5])

    sim.add_gate(Gates.Swap([1, 2], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0, 0, 0.5, 0, 0, 0.5])

    sim.add_gate(Gates.H([2], n_qbits))
    assert_arrays(sim.print_statevector(), [Gates.SQRT_H, Gates.SQRT_E, 0, Gates.SQRT_E, 0, Gates.SQRT_E, 0, -Gates.SQRT_E])

    sim.add_gate(Gates.Z([0], n_qbits))
    sim.add_gate(Gates.S([1], n_qbits))
    sim.add_gate(Gates.T([2], n_qbits))
    assert_arrays(sim.print_statevector(), [Gates.SQRT_H, -Gates.SQRT_E, 0, -Gates.SQRT_E * 1j, 0, -0.25 - 0.25j, 0, -0.25 + 0.25j])

    sim.print_sim()

def test_bell_state():
    n_qbits = 2
    sim = Sim(n_qbits, "Test: Bell State", True)
    print("")

    sim.add_gate(Gates.H([0], n_qbits))
    assert_arrays(sim.print_statevector(), [Gates.SQRT_H, Gates.SQRT_H, 0, 0])

    sim.add_gate(Gates.CX([0, 1], n_qbits))
    assert_arrays(sim.print_statevector(), [Gates.SQRT_H, 0, 0, Gates.SQRT_H])

    sim.print_sim()

def test_hadamard():
    n_qbits = 2
    sim = Sim(n_qbits, "Test: Hadamard", True)
    print("")

    sim.add_gate(Gates.H([0], n_qbits))
    assert_arrays(sim.print_statevector(), [Gates.SQRT_H, Gates.SQRT_H, 0, 0])

    sim.add_gate(Gates.H([1], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0.5])

    sim.print_sim()

def test_3():
    n_qbits = 3
    sim = Sim(n_qbits, "Test 3", True)
    print("")

    sim.add_gate(Gates.H([0], n_qbits))
    sim.add_gate(Gates.H([1], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0])

    sim.add_gate(Gates.CX([1, 2], n_qbits))
    assert_arrays(sim.print_statevector(), [0.5, 0.5, 0, 0, 0, 0, 0.5, 0.5])

    sim.print_sim()

def assert_arrays(arr1, arr2):
    assert arr1 == pytest.approx(arr2, margin)

def assert_matrices(a, b):
    np.testing.assert_array_almost_equal_nulp(a, b)

def main():
    test_statevector()
    test_bell_state()
    test_hadamard()
    test_3()


if __name__ == "__main__":
    np.set_printoptions(threshold=sys.maxsize)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
    main()