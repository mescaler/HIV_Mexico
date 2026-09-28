import random
from collections import defaultdict
from Bio import SeqIO

# Set random seed for reproducibility
random.seed(42)

# Input and output files
input_fasta = "input.fasta"
output_fasta = "R1.fasta"

# Desired subsample counts: {year: {location: count}}
subsample_counts = {
    '2019': {'ST_1': 310, 'ST_2': 229},
    '2020': {'ST_1': 166, 'ST_2': 111},
    '2021': {'ST_1': 288, 'ST_2': 214},
    '2022': {'ST_1': 314, 'ST_2': 264}
}

# Step 1: Parse and group sequences by year and location
seq_dict = defaultdict(list)

for record in SeqIO.parse(input_fasta, "fasta"):
    header_parts = record.description.split("|")
    if len(header_parts) < 5:
        continue  # Skip malformed headers

    location = header_parts[2]
    date = header_parts[4]
    year = date.split("/")[0]

    if year in subsample_counts and location in subsample_counts[year]:
        seq_dict[(year, location)].append(record)

# Step 2: Randomly sample sequences according to target counts
final_sample = []

for year in subsample_counts:
    for location in subsample_counts[year]:
        available = seq_dict.get((year, location), [])
        n_desired = subsample_counts[year][location]

        if len(available) < n_desired:
            raise ValueError(f"Not enough sequences for {location} in {year}: requested {n_desired}, available {len(available)}")

        sampled = random.sample(available, n_desired)
        final_sample.extend(sampled)

# Step 3: Write the subsampled sequences to output file
SeqIO.write(final_sample, output_fasta, "fasta")

print(f"Subsampled {len(final_sample)} sequences to {output_fasta}")