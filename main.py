import numpy as np
import sys
from gates import Gates
from sim import Sim

def sim_demo():
    n_qbits = 3
    sim = Sim(n_qbits, "Demo")

    sim.add_gate(Gates.H([0], n_qbits))
    sim.add_gate(Gates.H([1], n_qbits))
    #sim.add_gate(Gates.GATE_H([2], n_qbits))

    sim.add_gate(Gates.Toffoli([0, 1, 2], n_qbits))

    sim.print_sim()
    sim.print_prob_dist()
    sim.print_statevector()


# Bell State Simulation
def sim_bell_state():
    n_qbits = 2
    sim = Sim(n_qbits, "Bell State")

    sim.add_gate(Gates.H([0], n_qbits))
    sim.add_gate(Gates.CX([0, 1], n_qbits))

    sim.print_sim()
    sim.print_prob_dist()


def main():
    #sim_bell_state()
    sim_demo()


if __name__ == "__main__":
    np.set_printoptions(threshold=sys.maxsize)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
    main()