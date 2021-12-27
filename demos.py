import numpy as np
import sys
from gate import Gate
from sim import Sim


def demo_basic():
    sim = Sim(3, "Demo: Basic")

    sim.add_gate(Gate.H(), [0])
    sim.add_gate(Gate.H(), [1])
    sim.add_gate(Gate.Toffoli(), [0, 1, 2])
    sim.add_gate(Gate.Swap(), [1, 2])
    sim.add_gate(Gate.H(), [2])

    sim.print_sim()
    sim.print_prob_dist()
    sim.print_statevector()


def demo_control():
    sim = Sim(5, "Demo: Control")

    # Multi-controlled gate
    gate_target = Gate.H()
    num_control_bits = 3
    gate_control = Gate.MC(gate_target, num_control_bits)
    sim.add_gate(gate_control, [0, 3, 4, 1])

    sim.add_gate(Gate.CSwap(), [1, 0, 4])

    sim.print_sim()


def demo_custom_gate():
    sim = Sim(5, "Demo: Custom gates")

    sim_to_be_gate = Sim(2)
    sim_to_be_gate.add_gate(Gate.H(), [0])
    sim_to_be_gate.add_gate(Gate.CX(), [0, 1])

    custom_gate = sim_to_be_gate.to_gate("G1")
    print(f"Expected # of qubits for custom gate: {custom_gate.expected_qbits}")

    sim.add_gate(custom_gate, [0, 4])
    sim.add_gate(custom_gate, [2, 1])

    sim.print_sim()


def demo_gates_with_args():
    sim = Sim(5, "Demo: Gates with arguments")

    theta_045 = 1 * np.pi / 4
    theta_135 = 3 * np.pi / 4
    theta_225 = 5 * np.pi / 4
    theta_315 = 7 * np.pi / 4

    sim.add_gate(Gate.P(theta_045), [0])
    sim.add_gate(Gate.Rz(theta_135), [1])
    sim.add_gate(Gate.Ryy(theta_225), [2, 0])
    sim.add_gate(Gate.CP(theta_315), [3, 4])

    sim.print_sim()


def demo_multi_targeted_gates():
    sim = Sim(3, "Demo: Multi-targeted gates")

    sim.add_gate(Gate.MT(Gate.H(), 3), [0, 1, 2])
    sim.add_gate(Gate.MT(Gate.T(), 3), [0, 1, 2])
    sim.add_gate(Gate.MT(Gate.H(), 3), [0, 1, 2])
    sim.add_gate(Gate.MT(Gate.T(), 3), [0, 1, 2])

    sim.print_sim()
    sim.print_prob_dist()
    sim.print_statevector()


# Bell State Simulation
def demo_bell_state():
    sim = Sim(2, "Bell State")

    sim.add_gate(Gate.H(), [0])
    sim.add_gate(Gate.CX(), [0, 1])

    sim.print_sim()


# Grover's search algorithm
def demo_grover():
    n_qbits = 3
    i_qbits = [i for i in range(n_qbits)]
    sim = Sim(n_qbits, "Grover's Search Algorithm")

    # Equal superposition
    sim.add_gate(Gate.MT(Gate.H(), n_qbits), i_qbits)

    # Oracle for |101> and |110>
    sim.add_gate(Gate.MC(Gate.Z(), 1), [0, 1])
    sim.add_gate(Gate.MC(Gate.Z(), 1), [0, 2])

    # Diffuser
    gate_diffuser = general_diffuser(n_qbits)
    sim.add_gate(gate_diffuser, i_qbits)

    sim.print_sim()
    sim.print_prob_dist()


def general_diffuser(n_qbits):
    sim = Sim(n_qbits)

    i_qbits = [i for i in range(n_qbits)]
    sim.add_gate(Gate.MT(Gate.H(), n_qbits), i_qbits)
    sim.add_gate(Gate.MT(Gate.X(), n_qbits), i_qbits)
    sim.add_gate(Gate.MC(Gate.Z(), n_qbits - 1), i_qbits)
    sim.add_gate(Gate.MT(Gate.X(), n_qbits), i_qbits)
    sim.add_gate(Gate.MT(Gate.H(), n_qbits), i_qbits)

    sim.print_sim()
    return sim.to_gate("Dif")


def main():
    # demo_bell_state()
    # demo_basic()
    # demo_control()
    # demo_custom_gate()
    # demo_gates_with_args()
    # demo_multi_targeted_gates()
    demo_grover()


if __name__ == "__main__":
    np.set_printoptions(threshold=sys.maxsize)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
    main()
