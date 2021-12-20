import numpy as np
import sys
from gates import Gates
from sim import Sim

# Bell State Simulation
def sim_bell_state():
    sim = Sim(2, "Bell State")

    sim.add_gate([0], Gates.H)
    sim.add_gate([0, 1], Gates.CX)

    sim.print_sim()
    sim.print_prob_dist()


"""
    Representing a 4-bit CCC-Not gate with Toffoli gates
    CCC-Not gate: Control: [0, 1, 2] -> Controlled: [3]
    i_qbit [4] is an auxiliary qubit
"""


def sim_cccnot_with_toffoli():
    sim = Sim(5, "CCC-Not gate with Toffoli gates")

    sim.add_gate([0, 1, 4], Gates.TOFFOLI)
    sim.add_gate([2, 4, 3], Gates.TOFFOLI)

    sim.print_sim()

def sim_debug():
    sim = Sim(3, "Debug")
    print("")

    #sim.print_prob_dist()
    #sim.print_statevector()

    #sim.print_sim()


def main():
    #sim_bell_state()
    #sim_cccnot_with_toffoli()
    sim_debug()


if __name__ == "__main__":
    np.set_printoptions(threshold=sys.maxsize)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
    main()