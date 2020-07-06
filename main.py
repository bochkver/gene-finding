#!/usr/bin/python
import re
from fractions import Fraction
from enum import Enum, auto
from collections import defaultdict

class States(Enum):
    start = auto()
    non_coding_a = auto()
    non_coding_c = auto()
    non_coding_g = auto()
    non_coding_t = auto()
    start_codon_a = auto()
    start_codon_c = auto()
    start_codon_g = auto()
    start_codon_t = auto()
    start_codon_2 = auto()
    start_codon_3 = auto()
    internal_codon_1_a = auto()
    internal_codon_1_c = auto()
    internal_codon_1_g = auto()
    internal_codon_1_t = auto()
    internal_codon_2_a = auto()
    internal_codon_2_c = auto()
    internal_codon_2_g = auto()
    internal_codon_2_t = auto()
    internal_codon_3_a = auto()
    internal_codon_3_c = auto()
    internal_codon_3_g = auto()
    internal_codon_3_t = auto()
    stop_codon_1_t = auto()
    stop_codon_2_a = auto()
    stop_codon_2_g = auto()
    stop_codon_3_a = auto()
    stop_codon_3_g = auto()
    end = auto()


def non_coding_state_from_symbol(s):
    if s == "a": return States.non_coding_a
    if s == "c": return States.non_coding_c
    if s == "g": return States.non_coding_g
    if s == "t": return States.non_coding_t


def start_codon_from_gen(gen):
    if gen == "atg":
        return [States.start_codon_a, States.start_codon_2, States.start_codon_3]
    if gen == "ctg":
        return [States.start_codon_c, States.start_codon_2, States.start_codon_3]
    if gen == "gtg":
        return [States.start_codon_g, States.start_codon_2, States.start_codon_3]
    if gen == "ttg":
        return [States.start_codon_t, States.start_codon_2, States.start_codon_3]


def internal_codon_from_gen(gen):
    state = []
    if gen[0] == "a": state.append(States.internal_codon_1_a)
    if gen[0] == "c": state.append(States.internal_codon_1_c)
    if gen[0] == "g": state.append(States.internal_codon_1_g)
    if gen[0] == "t": state.append(States.internal_codon_1_t)
    if gen[1] == "a": state.append(States.internal_codon_2_a)
    if gen[1] == "c": state.append(States.internal_codon_2_c)
    if gen[1] == "g": state.append(States.internal_codon_2_g)
    if gen[1] == "t": state.append(States.internal_codon_2_t)
    if gen[2] == "a": state.append(States.internal_codon_3_a)
    if gen[2] == "c": state.append(States.internal_codon_3_c)
    if gen[2] == "g": state.append(States.internal_codon_3_g)
    if gen[2] == "t": state.append(States.internal_codon_3_t)
    return state


def stop_codon_from_gen(gen):
    if gen == "taa":
        return [States.stop_codon_1_t, States.stop_codon_2_a, States.stop_codon_3_a]
    if gen == "tag":
        return [States.stop_codon_1_t, States.stop_codon_2_a, States.stop_codon_3_g]
    if gen == "tga":
        return [States.stop_codon_1_t, States.stop_codon_2_g, States.stop_codon_3_a]


def train(sequence, genes):
    # make states_list for sequence
    states_list = [non_coding_state_from_symbol(s) for s in sequence]
    for start, end in genes:
        states_list[start:start + 3] = start_codon_from_gen(sequence[start:start + 3])
        for i in range(start + 3, end - 2, 3):
            states_list[i:i + 3] = internal_codon_from_gen(sequence[i:i + 3])
        states_list[end - 2:end + 1] = stop_codon_from_gen(sequence[end - 2:end + 1])

    # fill states_transition matrix
    states_transition[States.start][states_list[0]] += 1
    states_transition[states_list[-1]][States.end] += 1
    for previous, current in zip(states_list, states_list[1:]):
        states_transition[previous][current] += 1

    # fill states_output
    for state, output in zip(states_list, sequence):
        states_output[state][output] += 1


def transition_probability(from_state, to_state):
    if to_state not in states_transition[from_state]:
        return Fraction(0, 1)
    return states_transition[from_state][to_state]


# by https://en.wikipedia.org/wiki/Viterbi_algorithm
def viterbi(sequence):
    V = [{}]
    for state in States:
        V[0][state] = {"prob": Fraction(0, 1), "prev": None}
    V[0][States.start] = {"prob": Fraction(1, 1), "prev": None}
    states = list(States)

    # Run Viterbi when t > 0
    sequence = "*" + sequence # shift sequence by 1 to right because of silent start state
    for t in range(1, len(sequence)):
        V.append({})
        for state in states:
            V[t][state] = {"prob": Fraction(0, 1), "prev": None}
        for prev_st, data in V[t - 1].items():
            if data["prob"] == 0:
                continue
            for state, probability in states_transition[prev_st].items():
                tr_prob = data["prob"] * probability * states_output[state][sequence[t]]
                state_data = V[t][state]
                if tr_prob > state_data["prob"]:
                    V[t][state] = {"prob": tr_prob, "prev": prev_st}
    # end state
    V.append({})
    V[-1][States.end] = {"prob": Fraction(0, 1), "prev": None}
    for prev_st, data in V[-2].items():
        if data["prob"] == 0:
            continue
        if States.end not in states_transition[prev_st]:
            continue
        tr_prob = data["prob"] * states_transition[prev_st][States.end]
        state_data = V[-1][States.end]
        if tr_prob > state_data["prob"]:
            V[-1][States.end] = {"prob": tr_prob, "prev": prev_st}
    opt = [States.end]
    previous = States.end
    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    return opt


ecoli_file = open("sequenceOfEcoliStrainM54.txt", "r")
training_file = open("train.seqs", "r")
testing_file = open("test.seqs", "r")
ecoli = ""
for line in ecoli_file:
    ecoli += line[:-1]

# states is list of dictionaries, where ech element represent how often i-st state transit to other states
# dictionary key is destination state and value is number of transits
states_transition = {state: defaultdict(int) for state in States}
states_output = {state: {"a": 0, "c": 0, "g": 0, "t": 0} for state in States}

for line in training_file:
    indexes = re.split(" |\t|,|\[|\]", line[:-1])
    indexes = filter(lambda x: len(x), indexes)
    indexes = list(map(lambda x: int(x) - 1, indexes))
    start = indexes[0]
    end = indexes[1]
    # correct gene indexes train sequence
    indexes = list(map(lambda x: x - start, indexes))
    genes = zip(indexes[2::2], indexes[3::2])
    train(ecoli[start:end], genes)

# normalize states_trans_matrix , states_output
for state in States:
    norm_factor = sum(states_transition[state].values())
    states_transition[state] = {key: Fraction(value, norm_factor) for key, value in states_transition[state].items()}
    norm_factor2 = sum(states_output[state].values())
    if norm_factor2 == 0:
        norm_factor2 = 1
    states_output[state] = {key: Fraction(value, norm_factor2) for key, value in states_output[state].items()}

# print in predictions file
predictions_file = open("predictions.txt", "w")
predicted_genes = 0
predicted_non_empty = 0
length_of_genes = 0
length_of_overlap = 0
length_of_genes_ground = 0
for line in testing_file:
    indexes = re.split(" |\t|,|\[|\]", line[:-1])
    indexes = filter(lambda x: len(x), indexes)
    indexes = list(map(lambda x: int(x) - 1, indexes))
    start = indexes[0]
    end = indexes[1]
    # correct gene indexes test sequence
    genes = zip(indexes[2::2], indexes[3::2])
    starts_ground = indexes[2::2]
    ends_ground = indexes[3::2]
    sequence = ecoli[start:end + 1]
    path = viterbi(sequence)
    starts = []
    ends = []
    # find all start and end of codon
    for idx, state in enumerate(path[1:-1]):
        if (state == States.start_codon_a
                or state == States.start_codon_c
                or state == States.start_codon_g
                or state == States.start_codon_t):
            starts.append(idx + start + 1)
        if (state == States.stop_codon_3_a or state == States.stop_codon_3_g):
            ends.append(idx + start + 1)
    output = [str(start + 1), str(end + 1)]
    genes = " ".join(["[" + str(s) + ", " + str(e) + "]" for s, e in zip(starts, ends)])
    output.append(genes)
    output_line = "\t".join(output) + "\n"
    predictions_file.write(output_line)
    predicted_genes += len(starts)

    for start_ground, end_ground in zip(starts_ground, ends_ground):
        length_of_genes_ground += end_ground - start_ground + 1

    for start, end in zip(starts, ends):
        length_of_genes += end - start + 1
        for start_ground, end_ground in zip(starts_ground, ends_ground):
            if (start <= start_ground and start_ground <= end):
                predicted_non_empty += 1
                length_of_overlap += end - start_ground + 1
                break
            if (start_ground <= start and start <= end_ground):
                predicted_non_empty += 1
                length_of_overlap += end_ground - start + 1
                break

# calculate accuracy
accuracy_file = open("accuracy.txt", "w")

# The number of predicted genes.
accuracy_file.write(str(predicted_genes))
accuracy_file.write("\n")

# The number of predicted genes that have nonempty overlap with at least one of the ground-truth genes.
accuracy_file.write(str(predicted_non_empty))
accuracy_file.write("\n")

# The total length of all genes that your model predicted.
accuracy_file.write(str(length_of_genes))
accuracy_file.write("\n")

# The total length of the overlap between predicted genes and ground-truth genes.
accuracy_file.write(str(length_of_overlap))
accuracy_file.write("\n")

# The Jaccard index of predicted genes and ground-truth genes.
# Consider each gene as a set of nucleotide indices in the E.coli strain.
# More on Jaccard index may be found on this Wikipedia page.

jaccard_index = length_of_overlap / (length_of_genes + length_of_genes_ground - length_of_overlap)
accuracy_file.write(str(jaccard_index))
