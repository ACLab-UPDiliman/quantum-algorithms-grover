# Author: Jeffrey A. Aborot
# Description: This quantum algorithm uses Grover's quantum algorithm as a scheme. It provides a quantum circuit design of the oracle for string matching.

# We follow the definition of Grover's quantum search algorithm in Explorations in Quantum Computing by William P. Collins, 2011 where Grover's operator is defined as Q=(-H)(I_s)(Hcross)(I_t).

# Library imports
from numpy import dot, kron, zeros, empty, random
from math import log, sqrt, floor, pi
from operators import X, Z, CNOT, CkNOT, W, H
import collections

#============================================================================================================
# Constants definition
ket_zero = [1,0]
ket_one = [0,1]

# Functions definition
def binary_to_vector(binary_state):
    if binary_state == 0:
        return [1,0]
    elif binary_state == 1:
        return [0,1]

def query_substrings(input_register_state):
    # This functions corresponds to the step of querying M-length substrings for T given the superposition of states corresponding to indices of T.

    return output_register_state


def shift_phase_of_t(index_amplitudes, substring_register_state, pattern_register_state):
    # This function corresponds to the application of U_f_prime operator to the superposition state
    # to mark the solution state t by shifting its phase.

    # Define scratch_register_state
    # scratch_register_state = [[[0,0],[0,1]],[[0,1],[1,0]],[[1,0],[1,1]],[[1,1],[0,0]]] =
    # substring_register_state = \
    #     [[[[1,0], [1,0]], [[1,0], [0,1]]],
    #      [[[1,0], [0,1]], [[0,1], [1,0]]],
    #      [[[0,1], [1,0]], [[0,1], [0,1]]],
    #      [[[0,1], [0,1]], [[1,0], [1,0]]]]

    # pattern_register_state = [[0, 1],
    #                           [1, 0]]
    # pattern_register_state = [[[1,0],[0,1]], [[0,1],[1,0]]]

    # scratch_register_state = \
    #     [[[[1,0], [1,0], [1,0]], [[1,0], [1,0], [1,0]]],
    #      [[[1,0], [1,0], [1,0]], [[1,0], [1,0], [1,0]]],
    #      [[[1,0], [1,0], [1,0]], [[1,0], [1,0], [1,0]]],
    #      [[[1,0], [1,0], [1,0]], [[1,0], [1,0], [1,0]]]]

    # Construct scratch register
    substrings_vectors = []
    for substring_index in xrange(N):
        substring_vectors = []
        for symbol_index in xrange(M):
            symbol_vectors = []
            for bit_index in xrange(symbol_bit_count+1):
                symbol_vectors.append(ket_zero)
            substring_vectors.append(symbol_vectors)
        substrings_vectors.append(substring_vectors)
    scratch_register_state = substrings_vectors


    # Mlog(|E|) CNOT operations: control: substring_register_state as control; target: scratch register (E is read as sigma denoted by Sigma)
    # print '>> Mlog(|E|) CNOT operations on substring and scratch register'

    for substring_index in xrange(0, N): # 0...N-1 (count of substrings in T)
        for symbol_index in xrange(0, M): # 0...M-1 (count of symbols in each substring)
            for bit_index in xrange(0, symbol_bit_count): # 0...log(|E|)-1 (count of bits in each symbol)

                # print 'substring_index: ', substring_index, ' symbol_index: ', symbol_index, ' bit_index: ', bit_index
                # print 'substring_register_state[', substring_index,'][', symbol_index,'][', bit_index,']: ', substring_register_state[substring_index][symbol_index][bit_index]
                # print 'scratch_register_state[', substring_index, '][', symbol_index,'][', bit_index, ']: ', scratch_register_state[substring_index][symbol_index][bit_index]

                state_vector = list(kron(substring_register_state[substring_index][symbol_index][bit_index], scratch_register_state[substring_index][symbol_index][bit_index]))
                new_state_vector = dot(CNOT(), state_vector)
                # print 'state_vector:', list(state_vector), ' new_state_vector:', list(new_state_vector)
                if (list(state_vector) != list(new_state_vector)):
                    # print 'state_vector ', state_vector, ' != new_state_vector ', new_state_vector
                    scratch_register_state[substring_index][symbol_index][bit_index] = list(dot(X(), scratch_register_state[substring_index][symbol_index][bit_index]))
                # else:
                    # print 'state_vector ', state_vector, ' == new_state_vector ', new_state_vector

                # print 'scratch_register_state[', substring_index, '][', symbol_index, '][', bit_index, ']: ', \
                # scratch_register_state[substring_index][symbol_index][bit_index]
                # print '==============================='

    # print 'substring register state: '
    # for i in xrange(N):
    #     for k in xrange(M):
    #         print substring_register_state[i][k]
    #
    # print 'scratch register state: '
    # for i in xrange(N):
    #     for k in xrange(M):
    #         print scratch_register_state[i][k],
    #     print

    # Mlog(|E|) CNOT operations: control: pattern_register_state as control; target: scratch register (E is read as sigma denoted by Sigma)
    # print '>> Mlog(|E|) CNOT operations on pattern and scratch register'
    for i in xrange(0, N): # 0...N-1 (count of substrings in T)
        for k in xrange(0, M): # 0...M-1 (count of symbols in each substring)
            for j in xrange(0, symbol_bit_count): # 0...log(|E|)-1 (count of bits in each symbol)
                substring_index = i
                symbol_index = k
                bit_index = j

                # print 'substring_index:', substring_index,' symbol_index: ', symbol_index, ' bit_index: ', bit_index
                # print 'pattern_register_state[', symbol_index,'][', bit_index,']: ', pattern_register_state[symbol_index][bit_index]
                # print 'scratch_register_state[', substring_index, '][', symbol_index,'][', bit_index, ']: ', scratch_register_state[substring_index][symbol_index][bit_index]

                state_vector = list(kron(pattern_register_state[symbol_index][bit_index], scratch_register_state[substring_index][symbol_index][bit_index]))
                new_state_vector = list(dot(CNOT(), state_vector))
                # print 'state_vector:', list(state_vector), ' new_state_vector:', list(new_state_vector)
                if (list(state_vector) != list(new_state_vector)):
                    # print 'state_vector ', state_vector, ' != new_state_vector ', new_state_vector
                    scratch_register_state[substring_index][symbol_index][bit_index] = list(dot(X(), scratch_register_state[substring_index][symbol_index][bit_index]))
                # else:
                    # print 'state_vector ', state_vector, ' == new_state_vector ', new_state_vector

                # print 'scratch_register_state[', substring_index, '][', symbol_index, '][', bit_index, ']: ', \
                # scratch_register_state[substring_index][symbol_index][bit_index]
                # print '==============================='

    # print 'pattern register state: '
    # for k in xrange(M):
    #     print pattern_register_state[k]
    #
    # print 'scratch register state: '
    # for i in xrange(N):
    #     for k in xrange(M):
    #         print scratch_register_state[i][k],
    #     print

    # M C^{log(|E|)}NOT operations: control: [k(log(|E|) + 1) + log(|E|)]-th qubits of scratch register; [k(log(|E|) + 1) + log(|E|)]-th qubits of scratch register
    # print '>> M C^{log(|E|)}NOT operations on scratch register'
    for substring_index in xrange(N):
        for symbol_index in xrange(M):
            scratch_register_state[substring_index][symbol_index].reverse()

            state_vector = scratch_register_state[substring_index][symbol_index][0] # add state of target bit to state vector
            for bit_index in xrange(1, symbol_bit_count+1): # construct the state vector for the control bits starting from bit nearest to the target bit
                scratch_register_state[substring_index][symbol_index][bit_index] = list(dot(X(),scratch_register_state[substring_index][symbol_index][bit_index])) # apply X to control bit
                state_vector = kron(scratch_register_state[substring_index][symbol_index][bit_index], state_vector) # add state of control bit to state vector
                scratch_register_state[substring_index][symbol_index][bit_index] = list(dot(X(), scratch_register_state[
                    substring_index][symbol_index][bit_index])) # apply X to control bit

            new_state_vector = dot(CkNOT(symbol_bit_count), state_vector)

            scratch_register_state[substring_index][symbol_index].reverse()

            if list(state_vector) != list(new_state_vector): # if all control bits are set to 1, apply X to target bit
                # print 'state_vector: ', state_vector, ' != new_state_vector: ', new_state_vector
                scratch_register_state[substring_index][symbol_index][symbol_bit_count] = list(dot(X(),scratch_register_state[substring_index][symbol_index][symbol_bit_count]))
            # else:
                # print 'state_vector: ', state_vector, ' == new_state_vector: ', new_state_vector

    # print 'scratch register state: '
    # for substring_index in xrange(N):
    #     for symbol_index in xrange(M):
    #         print scratch_register_state[substring_index][symbol_index],
    #     print

    # 1 C^{M}NOT operation: control: [k(log(|E|) + 1) + log(|E|)]-th qubits of scratch register; target: only qubit of output register
    # print '>> 1 C^{M}NOT operations on scratch register and output register'

    # Construct output register state
    # output_register_state = [[1,0],[1,0],[1,0],[1,0]]
    output_register_state = []
    for i in xrange(N):
        output_register_state.append(ket_zero)

    for substring_index in xrange(N):
        state_vector = scratch_register_state[substring_index][0][symbol_bit_count]
        for symbol_index in xrange(M):
            state_vector = kron(scratch_register_state[substring_index][symbol_index][symbol_bit_count], state_vector)

        new_state_vector = dot(CkNOT(M), state_vector)
        if list(state_vector) != list(new_state_vector):
            # print 'state_vector: ', state_vector, ' != new_state_vector: ', new_state_vector
            output_register_state[substring_index] = list(dot(X(), output_register_state[substring_index]))
        # else:
            # print 'state_vector: ', state_vector, ' == new_state_vector: ', new_state_vector

    # print 'output_register_state: ', output_register_state

    # 1 Z operation: only qubit of output register
    # print '>> 1 Z operation on output register'
    for substring_index in xrange(N):
        output_register_state[substring_index] = list(dot(Z(), output_register_state[substring_index]))
    # print 'output_register_state: ', output_register_state

    # index_amplitudes = [1.0 / 2, 1.0 / 2, 1.0 / 2, 1.0 / 2]
    # print 'index_amplitudes: ', index_amplitudes
    for substring_index in xrange(N):
        index_amplitudes[substring_index] = index_amplitudes[substring_index] * sum(output_register_state[substring_index])
    # print 'index_amplitudes: ', index_amplitudes
    return index_amplitudes

def shift_phase_of_0(state):
    # This function corresponds to the application of operator I_s = I - 2\vert 0 \rangle\langle 0 \vert into the index register.

    return state

def shift_phase_of_superposition(index_amplitudes):
    # print '>> (2|Psi><Psi| - I) operation on index register'
    index_amplitudes = dot(W(index_bit_count), index_amplitudes)
    # print 'index_amplitudes: ', index_amplitudes
    return index_amplitudes

def measure_state_of_index_register(amplitudes):
    # This function corresponds to the measurement of the state of the index register.
    # Perform probabilistic measurement on array of amplitudes and return output index.
    index = 0
    return index

def to_superposition(index_register_state):
    # print '>> H operation on index register'
    # for state in index_register_state:
    #     state.reverse()
    #     state_vector = binary_to_vector(state[0])
    #     for bit_index in xrange(len(state)-1):
    #         state_vector = kron(binary_to_vector(state[bit_index]), state_vector)
    #
    #     matrix = H()
    #     for i in xrange(index_bit_count-1):
    #         matrix = kron(H(),matrix)
    #
    #     print 'state_vector: ', state_vector
    #     print 'matrix:', matrix
    #     state_vector = dot(matrix,state_vector)
    #     print 'state_vector: ', state_vector

    state_vector = index_register_state[0][0]
    for i in xrange(1, index_bit_count) :
        state_vector = kron(index_register_state[0][i], state_vector)

    matrix = H()
    for i in xrange(index_bit_count - 1):
        matrix = kron(H(), matrix)

    # print 'state_vector: ', state_vector
    # print 'matrix:', matrix
    index_amplitudes = dot(matrix, state_vector)
    # print 'index_amplitudes: ', index_amplitudes

    return index_amplitudes

def initialize_index_register(index_register_state, index_bit_count):
    state_vectors = []
    for i in xrange(index_bit_count):
        state_vectors.append(ket_zero)
    return state_vectors

def encode_pattern_vectors(list_P, M, symbol_bit_count):
    symbols_vectors = []
    for symbol in list_P:
        symbol_vectors = []
        for bit in symbol_encoding[symbol]: # symbol_encoding returns a list [{0|1},{0|1}]
            symbol_vectors.append(binary_to_vector(bit))
        symbols_vectors.append(symbol_vectors)
    return symbols_vectors

def symbol_to_binary_list(symbol):
    return symbol_encoding[symbol]

def query_substrings_vectors(T, N, M):
    substrings_vectors = []
    for i in xrange(N): # for each substring in T
        substring = list_T[i:i+M]
        # print 'list_T[',i,':',i+M,']: ', substring
        substring_vectors = []
        for j in xrange(M): # for each symbol in each substring
            symbol_binary_encoding = symbol_to_binary_list(substring[j])
            symbol_vectors = []
            for bit in symbol_binary_encoding:
                symbol_vectors.append(binary_to_vector(bit))
            substring_vectors.append(symbol_vectors)
        substrings_vectors.append(substring_vectors)
    # print 'substring_vectors: ', substrings_vectors
    return substrings_vectors
#============================================================================================================
# Body
# NOTE: Encode states in vector format notation.
# NOTE: Perform computation as matrix operation on vector.

# Binary indexing of alphabet symbols: a=00, c=01, t=10, g=11
alphabet = ['a','c','t','g']
# symbol_encoding = {'a':[0,0],'c':[1,0],'t':[0,1],'g':[1,1], '-':[0,0]}
symbol_encoding = {'a':[0,0,0],'c':[1,0,0],'t':[0,1,0],'g':[1,1,0], '-':[0,0,1]}
alphabet_size = len(alphabet)
symbol_bit_count = int(log(alphabet_size, 2))+1
T = 'actggggggggggggg'
list_T = list(T)
N = len(T)
index_bit_count = int(log(N, 2))
P = 'actg'
list_P = list(P)
M = len(P)

for i in xrange(M-1):
    list_T.append('-')
# print 'list_T: ', list_T

# TODO: Add 1 qubit to accommodate filler symbol -.
# Initialize index register state
index_register_state = []
index_register_state.append(initialize_index_register(index_register_state, index_bit_count))
# print 'index_register_state: ', index_register_state

# Construct array of amplitudes of states in superposition
index_amplitudes = to_superposition(index_register_state)

# TODO: Implement simulation of a qRAM.
# Query substrings from QRAM and encode them into substring register as a superposition of states
substring_register_state = query_substrings_vectors(T, N, M)

# pattern_register_state = [[1,0],[0,1]] #= [[[1,0],[0,1]],[[0,1],[1,0]]] in vector form
pattern_register_state = encode_pattern_vectors(list_P, M, symbol_bit_count)
# print 'pattern_register_state: ', pattern_register_state

iter_count = int(floor((pi/4) * sqrt(N)))
# print 'iter_count: ', iter_count

for iter in xrange(iter_count):
    print '\n>> Grover iteration #', (iter + 1)

    print '\n>>>> Marking solution states: '
    index_amplitudes = shift_phase_of_t(index_amplitudes, substring_register_state, pattern_register_state)
    probabilities = [pow(i, 2) for i in index_amplitudes]
    print '>>>>>> amplitudes: ', index_amplitudes
    print '>>>>>> probabilities: ', probabilities
    print '>>>>>> total probabilities: ', sum(probabilities)

    print '\n>>>> Amplifying amplitude of solution states '
    # index_register_state = apply_H_adjoint(index_register_state)
    # index_register_state = shift_phase_of_0(index_register_state)
    # index_register_state = apply_H_negative(index_register_state)
    index_amplitudes = shift_phase_of_superposition(index_amplitudes)
    probabilities = [pow(i, 2) for i in index_amplitudes]
    print '>>>>>> amplitudes: ', index_amplitudes
    print '>>>>>> probabilities: ', probabilities
    print '>>>>>> total probabilities: ', sum(probabilities)

print '\n>> Measuring state of index register'
# Single sampling
print 'Output: {:d}'.format(random.choice(a=N,p=probabilities))

# Multiple sampling
# print 'i | count'
# index = measure_state_of_index_register(index_amplitudes)
# for i in xrange(N):
#     if i < 10:
#         print '0{:d}'.format(i), '| ', list(random.choice(N,100,True,probabilities)).count(i)
#     else:
#         print '{:d}'.format(i), '| ', list(random.choice(N, 100, True, probabilities)).count(i)