import numpy as np
import sys
import pytest
from gates import Gates
from sim import Sim

margin = 0.000000001

def test_0():
    sim = Sim(3, "Test 0", True)
    print("")

    sim.add_gate([0], Gates.H)
    sim.add_gate([1], Gates.H)
    check_arrays(sim.print_statevector(),   [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0])

    sim.add_gate([0, 1, 2], Gates.TOFFOLI)
    check_arrays(sim.print_statevector(),   [0.5, 0.5, 0.5, 0, 0, 0, 0, 0.5])

    sim.add_gate([1, 2], Gates.SWAP)
    check_arrays(sim.print_statevector(),   [0.5, 0.5, 0, 0, 0.5, 0, 0, 0.5])

    sim.add_gate([2], Gates.H)
    check_arrays(sim.print_statevector(),   [Gates.SQRT_H, Gates.SQRT_E, 0, Gates.SQRT_E, 0, Gates.SQRT_E, 0, -Gates.SQRT_E])

    sim.add_gate([0], Gates.Z)
    sim.add_gate([1], Gates.S)
    sim.add_gate([2], Gates.T)
    check_arrays(sim.print_statevector(),   [Gates.SQRT_H, -Gates.SQRT_E, 0, -Gates.SQRT_E * 1j, 0, -0.25 - 0.25j, 0, -0.25 + 0.25j])

    sim.print_sim()

def test_1():
    sim = Sim(2, "Test 1", True)
    print("")

    sim.add_gate([0], Gates.H)
    check_arrays(sim.print_statevector(),   [Gates.SQRT_H, Gates.SQRT_H, 0, 0])

    sim.add_gate([0, 1], Gates.CX)
    check_arrays(sim.print_statevector(),   [Gates.SQRT_H, 0, 0, Gates.SQRT_H])

    sim.print_sim()

def test_hadamard():
    sim = Sim(2, "Test Hadamard", True)
    print("")

    sim.add_gate([0], Gates.H)
    check_arrays(sim.print_statevector(),   [Gates.SQRT_H, Gates.SQRT_H, 0, 0])

    sim.add_gate([1], Gates.H)
    check_arrays(sim.print_statevector(),   [0.5, 0.5, 0.5, 0.5])

    sim.print_sim()

def test_3():
    sim = Sim(3, "Test 3", True)
    print("")

    sim.add_gate([0], Gates.H)
    sim.add_gate([1], Gates.H)
    check_arrays(sim.print_statevector(),   [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0])

    sim.add_gate([1, 2], Gates.CX)
    check_arrays(sim.print_statevector(),   [0.5, 0.5, 0, 0, 0, 0, 0.5, 0.5])

    sim.print_sim()

def test_gates_args():
    sim = Sim(3, "Test Gates with Arguments", True)
    print("")

    sim.print_sim()

def check_arrays(arr1, arr2):
    assert arr1 == pytest.approx(arr2, margin)

def main():
    test_0()
    test_1()
    test_hadamard()
    test_3()
    test_gates_args()


if __name__ == "__main__":
    np.set_printoptions(threshold=sys.maxsize)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
    main()