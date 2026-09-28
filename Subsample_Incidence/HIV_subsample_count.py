from Bio import SeqIO
from collections import defaultdict

fasta_path = "R1.fasta"
counts = defaultdict(lambda: defaultdict(int))

# Parse fasta and count
for record in SeqIO.parse(fasta_path, "fasta"):
    header_parts = record.description.split("|")
    location = header_parts[2]       # e.g. 'ST_1' or 'ST_2'
    year = header_parts[4].split("/")[0]  # e.g. '2019' from '2019/07/17'
    counts[year][location] += 1

# Get all years and locations sorted
years = sorted(counts.keys())
locations = sorted({loc for year in years for loc in counts[year].keys()})

# Prepare totals
col_sums = {loc: 0 for loc in locations}
row_sums = {}

# Print header
print(f"{'year':<8}" + "".join(f"{loc:>6}" for loc in locations) + "  Sum")

# Print rows with counts and row sums
for year in years:
    row_total = 0
    row_values = []
    for loc in locations:
        val = counts[year].get(loc, 0)
        row_values.append(val)
        row_total += val
        col_sums[loc] += val
    row_sums[year] = row_total
    print(f"{year:<8}" + "".join(f"{v:>6}" for v in row_values) + f"  {row_total}")

# Print column sums row
total_all = sum(row_sums.values())
print(f"{'Sum':<8}" + "".join(f"{col_sums[loc]:>6}" for loc in locations) + f"  {total_all}\n")

print(f"Total sequences: {total_all}")