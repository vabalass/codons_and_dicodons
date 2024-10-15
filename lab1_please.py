import os
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import numpy as np

starts = ['ATG']
stops = ['TAA', 'TAG', 'TGA']
sequences_dir = "data"

def extract_frames(sequence_data):
    all_frames = []
    rev_comp_sequence = sequence_data.reverse_complement()
    for offset in range(3):
        all_frames.append(sequence_data[offset:])
        all_frames.append(rev_comp_sequence[offset:])
    return all_frames

def locate_coding_sequences(frame_sequence):
    coding_regions = []
    seq_length = len(frame_sequence)

    for index in range(0, seq_length - 2, 3):
        current_codon = frame_sequence[index:index + 3]
        if current_codon in starts:
            for j in range(index + 3, seq_length - 2, 3):
                potential_stop = frame_sequence[j:j + 3]
                if potential_stop in stops:
                    extracted_seq = frame_sequence[index:j + 3]
                    if len(extracted_seq) > 100:
                        coding_regions.append(extracted_seq)
                    break
    return coding_regions

def analyze_coding_regions(coding_sequences):
    codon_freq_map = defaultdict(int)
    dicodon_freq_map = defaultdict(int)
    codon_total = 0
    dicodon_total = 0

    for sequence in coding_sequences:
        for k in range(0, len(sequence) - 2, 3):
            codon = sequence[k:k + 3]
            protein = str(Seq(codon).translate())
            if protein != '*':
                codon_freq_map[protein] += 1
                codon_total += 1

            if k + 6 <= len(sequence):
                dicodon = sequence[k:k + 6]
                double_protein = str(Seq(dicodon).translate())
                if '*' not in double_protein:
                    dicodon_freq_map[double_protein] += 1
                    dicodon_total += 1

    codon_frequencies = {k: v / codon_total for k, v in codon_freq_map.items()} if codon_total > 0 else {}
    dicodon_frequencies = {k: v / dicodon_total for k, v in dicodon_freq_map.items()} if dicodon_total > 0 else {}

    return codon_frequencies, dicodon_frequencies

def generate_distance_matrix(frequency_data):
    sequence_keys = list(frequency_data.keys())
    num_sequences = len(sequence_keys)
    unique_codons = set()

    for values in frequency_data.values():
        unique_codons.update(values.keys())

    distance_matrix = np.zeros((num_sequences, num_sequences))

    for x in range(num_sequences):
        for y in range(num_sequences):
            if x != y:
                codon_vec_x = np.array([frequency_data[sequence_keys[x]].get(codon, 0) for codon in unique_codons])
                codon_vec_y = np.array([frequency_data[sequence_keys[y]].get(codon, 0) for codon in unique_codons])
                squared_difference = np.abs(codon_vec_x ** 2 - codon_vec_y ** 2)
                distance_matrix[x][y] = np.sum(squared_difference)

    return distance_matrix

def generate_dicodon_matrix(dicodon_data):
    sequence_list = list(dicodon_data.keys())
    sequence_count = len(sequence_list)
    unique_dicodons = set()

    for frequencies in dicodon_data.values():
        unique_dicodons.update(frequencies.keys())

    matrix = np.zeros((sequence_count, sequence_count))

    for m in range(sequence_count):
        for n in range(sequence_count):
            if m != n:
                dicodon_vec_m = np.array([dicodon_data[sequence_list[m]].get(dicodon, 0) for dicodon in unique_dicodons])
                dicodon_vec_n = np.array([dicodon_data[sequence_list[n]].get(dicodon, 0) for dicodon in unique_dicodons])
                diff_squared = np.abs(dicodon_vec_m ** 2 - dicodon_vec_n ** 2)
                matrix[m][n] = np.sum(diff_squared)

    return matrix

def output_to_phylip(file_output_name, distance_matrix, sequence_names):
    with open(file_output_name, 'w') as phylip_file:
        phylip_file.write(f"{len(sequence_names)}\n")
        for i in range(len(sequence_names)):
            distance_string = " ".join(f"{distance_matrix[i][j]:.3f}" for j in range(len(sequence_names)))
            name = sequence_names[i][:10].ljust(10)
            phylip_file.write(f"{name} {distance_string}\n")

def retrieve_top_items(frequency_dictionary, limit=10):
    sorted_freq = sorted(frequency_dictionary.items(), key=lambda item: item[1], reverse=True)
    return sorted_freq[:limit]

def analyze_fasta_data():
    codon_data = defaultdict(lambda: defaultdict(int))
    dicodon_data = defaultdict(lambda: defaultdict(int))
    virus_codon_count_mammalian = defaultdict(int)
    virus_codon_count_bacterial = defaultdict(int)
    dicodon_count_mammalian = defaultdict(int)
    dicodon_count_bacterial = defaultdict(int)

    for fasta_file in os.listdir(sequences_dir):
        if fasta_file.endswith(".fasta"):
            file_location = os.path.join(sequences_dir, fasta_file)
            print(f"Analyzing file: {fasta_file}")

            for fasta_record in SeqIO.parse(file_location, "fasta"):
                nucleotide_seq = fasta_record.seq
                all_frames = extract_frames(nucleotide_seq)

                all_coding_regions = []
                for single_frame in all_frames:
                    coding_region_set = locate_coding_sequences(single_frame)
                    all_coding_regions.extend(coding_region_set)

                codon_freqs, dicodon_freqs = analyze_coding_regions(all_coding_regions)

                if "mamalian" in fasta_file.lower():
                    for codon_key, freq_value in codon_freqs.items():
                        virus_codon_count_mammalian[codon_key] += freq_value
                    for dicodon_key, freq_value in dicodon_freqs.items():
                        dicodon_count_mammalian[dicodon_key] += freq_value
                elif "bacterial" in fasta_file.lower():
                    for codon_key, freq_value in codon_freqs.items():
                        virus_codon_count_bacterial[codon_key] += freq_value
                    for dicodon_key, freq_value in dicodon_freqs.items():
                        dicodon_count_bacterial[dicodon_key] += freq_value

                for codon_key, freq_value in codon_freqs.items():
                    codon_data[fasta_file][codon_key] += freq_value
                for dicodon_key, freq_value in dicodon_freqs.items():
                    dicodon_data[fasta_file][dicodon_key] += freq_value

    codon_matrix = generate_distance_matrix(codon_data)
    dicodon_matrix = generate_dicodon_matrix(dicodon_data)

    output_to_phylip("codon_phylip_output.phy", codon_matrix, list(codon_data.keys()))
    output_to_phylip("dicodon_phylip_output.phy", dicodon_matrix, list(dicodon_data.keys()))

    print("\nTop 10 Codons in Mammalian Viruses:")
    for codon, count in retrieve_top_items(virus_codon_count_mammalian):
        print(f"{codon}: {count:.4f}")

    print("\nTop 10 Codons in Bacterial Viruses:")
    for codon, count in retrieve_top_items(virus_codon_count_bacterial):
        print(f"{codon}: {count:.4f}")

    print("\nTop 10 Dicodons in Mammalian Viruses:")
    for dicodon, count in retrieve_top_items(dicodon_count_mammalian):
        print(f"{dicodon}: {count:.4f}")

    print("\nTop 10 Dicodons in Bacterial Viruses:")
    for dicodon, count in retrieve_top_items(dicodon_count_bacterial):
        print(f"{dicodon}: {count:.4f}")

if __name__ == "__main__":
    analyze_fasta_data()
