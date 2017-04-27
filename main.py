# Author: Jeffrey A. Aborot
# Description: This quantum algorithm uses Grover's quantum algorithm as a scheme. It provides a quantum circuit design of the oracle for string matching.

# We follow the definition of Grover's quantum search algorithm in Explorations in Quantum Computing by William P. Collins, 2011 where Grover's operator is defined as Q=(-H)(I_s)(Hcross)(I_t).

#============================================================================================================
# Functions definition
def query_substrings(input_register_state):
    # This functions corresponds to the step of querying M-length substrings for T given the superposition of states corresponding to indices of T.

    return output_register_state


def shift_phase_of_t(index_register_state, substring_register_state, pattern_register_state):
    # This function corresponds to the application of U_f_prime operator to the superposition state
    # to mark the solution state t by shifting its phase.
    scratch_register_state # express state in vector format notation
    output_register_state # express state in vector format notation

    return index_state

def apply_H_adjoint(state):
    # This function applies the operator (H^{\cross})^{\otimes \log(N)} to the index register.

    return state

def shift_phase_of_0(state):
    # This function corresponds to the application of operator I_s = I - 2\vert 0 \rangle\langle 0 \vert into the index register.

    return state

def apply_H_negative(state):
    # This function applies the operator -(H)^{\otimes \log(N)} to the index register.

    return state

def measure_state_of_index_register(state):
    # This function corresponds to the measurement of the state of the index register.

    return state

#============================================================================================================
# Body
# NOTE: Encode states in vector format notation.
# NOTE: Perform computation as matrix operation on vector.

# Binary indexing of alphabet symbols: a=00, c=01, t=10, g=11
# T = actg = 00,01,10,11
# P = ct = 01,10
# substrings of T:
# ac = 00,01
# ct = 01,10
# tg = 10,11
# g- = 11,00
# TODO: Add 1 qubit to accommodate filler symbol -.
index_register_state = [[0,0]] #= [[1,0],[1,0]] in vector form
#substring_register_state = [[[0,0],[0,1]],[[0,1],[1,0]],[[1,0],[1,1]],[[1,1],[0,0]]]
pattern_register_state = [[0,1],[1,0]] #= 1/sqrt(2) * [[[1,0],[0,1]],[[0,1],[1,0]]] in vector form

# TODO: Implement superposition operation on 00 state.
index_register_state = [[0,0],[0,1],[1,0],[1,1]] #= 1/2 * [[[1,0],[1,0]],[[1,0],[0,1]],[[0,1].[1,0]],[[0,1].[0,1]]] in vector form

# TODO: Implement simulation of a qRAM.
substring_register_state = query_substrings(index_register_state)

# TODO: Implement marking operation on superposition state.
index_register_state = shift_phase_of_t(substring_register_state, pattern_register_state)

# TODO: Implement amplitude amplification of marked state.
index_register_state = apply_H_adjoint(index_register_state)
index_register_state = shift_phase_of_0(index_register_state)
index_register_state = apply_H_negative(index_register_state)

# TODO: Implement measurement operation.
t = measure_state_of_index_register(index_register_state)