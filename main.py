import numpy as np

SQRT_H = 1 / np.sqrt(2)

# Qubit base vectors
STATE_0 = np.array([1, 0])
STATE_1 = np.array([0, 1])

# Classical bit states
CSTATE_0 = 0
CSTATE_1 = 1

# Gate labels
GATE_I = "I"
GATE_X = "X"
GATE_Y = "Y"
GATE_Z = "Z"
GATE_H = "H"
GATE_S = "S"
GATE_T = "T"
GATE_CNOT = "C"
GATE_TOFFOLI = "O"
GATE_MEASUREMENT = "M"
GATE_NOTHING = "-"

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

GATE_CNOT_M = np.array(
    [[1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 0, 1],
     [0, 0, 1, 0]])

GATE_TOFFOLI_M = np.array(
    [[1, 0, 0, 0, 0, 0, 0, 0],
     [0, 1, 0, 0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0, 0, 0, 0],
     [0, 0, 0, 1, 0, 0, 0, 0],
     [0, 0, 0, 0, 1, 0, 0, 0],
     [0, 0, 0, 0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 1, 0]])

"""
    All qubits are initialized to the state [1, 0] at the start
    
    t: int
    i_qbit: int
    i_qbits: int array
    gate: gate label (GATE_*)
    state: int array
    num_samples: int
"""


class Sim:
    n_qbits = 0
    max_t = 0

    """
        states of qubits
        
        [i_qbit] = state
    """
    qbits = {}

    """
        list of gates (np arrays) of type "GATE_*"
        
        [t] = [i_qbits, gate]
    """
    gates = {}

    """
        initializes the simulator with n_qbits qubits
            and n_cbits classical bits
    """

    def __init__(self, n_qbits):
        self.n_qbits = n_qbits

        for i_qbit in range(self.n_qbits):
            self.qbits[i_qbit] = STATE_0

    """
        adds the gate to the simulator at time t
            with the qbits listed in i_qbits
    """

    def add_gate(self, i_qbits, gate):
        if self.__check_valid_i_qbits(i_qbits) != 0:
            print("Cannot add gate, invalid i_qbits: " + i_qbits.__repr__())
            return

        self.gates[self.max_t] = [i_qbits, gate]
        self.max_t += 1

    """
        adds a measurement for i_qbit as a "gate"
    """

    def add_measurement(self, i_qbit):
        if i_qbit >= self.n_qbits:
            print("Cannot add measurement, invalid i_qbit: " + i_qbit)

        self.gates[self.max_t] = [[i_qbit], GATE_MEASUREMENT]
        self.max_t += 1

    def __check_valid_i_qbits(self, i_qbits):
        for i_qbit in i_qbits:
            if i_qbit >= self.n_qbits:
                return -1

        return 0

    """
        run the simulator and display results
    """

    def run_sim(self, num_samples):
        if self.__validate_sim() > 0:
            return

        for i_sample in range(num_samples):
            self.__run_one_sample()

    """
        run the simulation once for a single sample
        save the results
    """

    def __run_one_sample(self):
        pass

    """
        ensure the simulator components are in a valid ordering
        - no gates or measurements added after a measurement
    """

    def __validate_sim(self):
        invalid_qbits = []

        for i_qbit in range(self.n_qbits):
            measurement_found = False

            for t in range(self.max_t):
                gate = self.gates[t]
                i_qbits = gate[0]
                gate_label = gate[1]

                if i_qbit in i_qbits:
                    if measurement_found:
                        invalid_qbits.append(i_qbit)
                        break

                    if gate_label is GATE_MEASUREMENT:
                        measurement_found = True

        self.__print_invalid_qbits(invalid_qbits)
        return len(invalid_qbits)

    def __print_invalid_qbits(self, invalid_qbits):
        if len(invalid_qbits) == 0:
            return

        print("Can't run simulation, invalid measurements at i_qbits: " + invalid_qbits.__repr__())

    """
        prints the quantum algorithm for this sim
        - doesn't identify inputs of 2-bit or 3-bit gates
    """

    def print_sim(self):
        for i_qbit in range(self.n_qbits):
            print("Q" + i_qbit.__repr__() + ": ", end="")
            for t in range(self.max_t):
                self.__print_single_qbit_at_t(t, i_qbit)
            print("")
        print("")

    def __print_single_qbit_at_t(self, t, i_qbit):
        gate = self.gates[t]
        i_qbits = gate[0]
        gate_label = gate[1]

        if i_qbit in i_qbits:
            print(gate_label, end="")
        else:
            print(GATE_NOTHING, end="")


"""
    Bell State Simulation
"""


def test_bell_state():
    print("Running test: Bell State")
    sim = Sim(2)

    sim.add_gate([0], GATE_H)
    sim.add_gate([0, 1], GATE_CNOT)

    sim.add_measurement(0)
    sim.add_measurement(1)

    sim.print_sim()
    sim.run_sim(100)


"""
    Example of an invalid simulator setup
    
    1. gates were added for qbits [0, 1, 2, 3, 4]
    2. measurements for qbits [0, 1, 3] were added
    3. gates were added for qbits [0, 1, 2]
    
    Since there were gates added for qbits [0, 1] after their
        measurements were taken, this simulator setup is invalid
"""


def test_invalid_measurement_placed():
    print("Running test: Invalid measurement placed")

    sim = Sim(5)

    sim.add_gate([0], GATE_S)
    sim.add_gate([1], GATE_X)
    sim.add_gate([2], GATE_Z)
    sim.add_gate([3, 4], GATE_CNOT)

    sim.add_measurement(0)
    sim.add_measurement(1)
    sim.add_measurement(3)

    sim.add_gate([0], GATE_Y)
    sim.add_gate([1], GATE_T)
    sim.add_gate([2], GATE_X)

    sim.print_sim()
    sim.run_sim(100)


def main():
    test_bell_state()
    test_invalid_measurement_placed()


if __name__ == "__main__":
    main()
