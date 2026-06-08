###########################################
# Genome Assembly as Shortest Superstring #
###########################################
from .overlap_graph import overlapZ

def assemble(sequences):
    sequence_names = list(sequences.keys())
    num_sequences = len(sequences)
    edges = {}

    min_length = len(sequences[sequence_names[0]]) / 2

    # Construct overlap graph
    for i in range(num_sequences):
        for j in range(num_sequences):
            if i != j:
                overlap_value = int(overlapZ(sequences[sequence_names[i]], sequences[sequence_names[j]]))
                if overlap_value >= min_length:
                    edges[(sequence_names[i], sequence_names[j])] = overlap_value

    contig_count = 0

    # Greedily remove edge corresponding to longest overlap
    while edges:
        (seq_name_i, seq_name_j), overlap_value = max(edges.items(), key=lambda x: x[1])

        seq_i = sequences[seq_name_i]
        seq_j = sequences[seq_name_j]

        contig_count += 1
        contig_name = ">contig_" + str(contig_count)
        sequences[contig_name] = seq_i + seq_j[overlap_value:]

        new_edges = {}
        for (a, b) in edges:
            if a == seq_name_i or a == seq_name_j or b == seq_name_i or b == seq_name_j:
                if b == seq_name_i and a != seq_name_j:
                    new_edges[(a, contig_name)] = edges[(a, b)]
                elif a == seq_name_j and b != seq_name_i:
                    new_edges[(contig_name, b)] = edges[(a, b)]
            else:
                new_edges[(a, b)] = edges[(a, b)]

        del sequences[seq_name_i]
        del sequences[seq_name_j]
        edges = new_edges

    if len(sequences) == 1:
        return sequences[list(sequences.keys())[0]] 
    else:
        return sequences