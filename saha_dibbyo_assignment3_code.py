import pandas as pd
import numpy as np

file_path = './gene_sequence.csv'
data = pd.read_csv(file_path)

def find_mismatches(sequence, target):
    mismatches = 0
    for i, j in zip(sequence, target):
        if i != j:
            mismatches +=1
    return mismatches

def shine_dalgarno(sequence):
    target = 'AGGAGG'
    minimum_mismatches = float('inf')
    best_position = None
    best_segment = None
    for i in range(17 - len(target) +1):
        segment = sequence[i:i+len(target)]
        mismatches = find_mismatches(segment, target)
        if mismatches < minimum_mismatches:
            minimum_mismatches = mismatches
            best_position = i
            best_segment = segment
        elif mismatches == minimum_mismatches:
            best_position = None
            best_segment = None
    return best_segment, minimum_mismatches, best_position

results = []
shine_dalgarno_sequences = []
for _, row in data.iterrows():
    gene = row['Gene']
    sequence = row['Sequence'].upper()
    shine_dalgarno_seq, mismatches, position = shine_dalgarno(sequence[:17])
    separation = position if position is not None else "None"
    if shine_dalgarno_seq is None:
        shine_dalgarno_seq = "None"
        mismatches = "None"
    results.append({'Gene': gene, 'Shine': shine_dalgarno_seq, 'Mismatches': mismatches, 'Separation':separation, '17 upstream bases': sequence[:17]})
    if shine_dalgarno_seq != "None":
        shine_dalgarno_sequences.append(shine_dalgarno_seq)
results_df = pd.DataFrame(results)
results_df.to_csv('saha_dibbyo_assignment3_spreadsheet.csv', index=False)

mismatches_values = [r['Mismatches'] for r in results if r['Mismatches'] != "None"]
separation_values = [r['Separation'] for r in results if r['Separation'] != "None"]

average_mismatches = np.mean(mismatches_values) if mismatches_values else None
print("Average Mismatches: ", average_mismatches)
standard_deviation_mismatches = np.std(mismatches_values) if mismatches_values else None
print("Standard Deviation Mismatches: ", standard_deviation_mismatches)
average_separation = np.mean(separation_values) if separation_values else None
print("Average Separation: ", average_separation)
standard_deviation_separation = np.std(separation_values) if separation_values else None
print("Standard Deviation Separation: ", standard_deviation_separation)

counts = [{} for _ in range(6)]
for seq in shine_dalgarno_sequences:
    for i, letter in enumerate(seq):
        if letter not in counts[i]:
            counts[i][letter] = 0
        counts[i][letter] += 1
consensus = ""
for count in counts:
    total = sum(count.values())
    for letter, count_temp in count.items():
        percentage = (count_temp / total) * 100
        print(f"Position: {counts.index(count) + 1}, Letter {letter}: {percentage:.2f}%")
    consensus_letter = max(count, key=count.get)
    consensus += consensus_letter
print("Consensus: ", consensus)