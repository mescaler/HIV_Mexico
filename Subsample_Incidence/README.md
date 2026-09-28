# Subsample Incidence HIV-GMC

Builds the epidemiologically-informed replica datasets used to test whether sampling bias between the two geographic locations (roughly 3:1 CDMX vs MEX) affects the results.

Each replica keeps all background sequences plus a randomly down-sampled fraction of the CIENI-CEC sequences (n = 1,896). The number kept per year per location is set by the lowest sequences-to-incidence ratio across the year × location cells, where incidence is new HIV cases per 100,000 inhabitants per year as officially reported by DGE, independently of CIENI-CEC. The per-cell targets are in Supplementary File 1, Supplementary Table 1.

**Associated with the preprint**

Proposed title: "HIV Transmission Dynamics Reveal Extensive Spatial Mixing Across Greater Mexico City"

Posted as:

Escalera-Zamudio M, López Ortiz E, García Morales C, Cruz-Bonilla E, Guerrero Flores S, Weaver S,
Matías Florentino M, Tapia Trejo D, Dávila Conn V, CIENI-CEC Consortium, Hernandez M, Dellicour S,
Vasylyeva TI, Wertheim JO, Ávila Ríos S.
*HIV Transmission Dynamics in Greater Mexico City are Shaped by Dense Spatial Mixing.*
medRxiv 2026.05.26.26354122; posted 27 May 2026.
<https://doi.org/10.64898/2026.05.26.26354122>

Contact: m.zamudio@ucl.ac.uk 

## Files

| File | What it is |
|------|------------|
| `HIV_subsample_incidence.py` | Reads `input.fasta`, keeps sequences from `ST_1`/`ST_2` in 2019–2022, randomly samples a fixed number per year × site, and writes `R1.fasta`. |
| `HIV_subsample_count.py` | Reads `R1.fasta` and prints a year × site count table with row and column totals. |
| `Supplementary_Table_1.csv` | Supplementary Table 1 with the SUBSAMPLE column recomputed as `round(incidence × 21.81)` |
| `NAME_DATE_REPS.xlsx` | Name/date tables for the replicas, one sheet each (`R1`–`R5`), 32,182 sequences per sheet (the 1,896 subsampled CIENI-CEC sequences plus all background sequences). Tip-date file for the time-scaled phylogeny. |

## Pipeline

```
input.fasta ──► HIV_subsample_incidence.py ──► R1.fasta ──► HIV_subsample_count.py ──► count table
                                                  │
                                                  └─► (combined with background sequences outside this folder)
                                                        ──► NAME_DATE_REPS.xlsx, sheet R1
```

### Expected FASTA header format

Both scripts split the header on `|` and expect:

```
<ID>|<country>|<site>|<cluster>|<YYYY/MM/DD>
  0      1        2       3          4
```

For example, `SEQID|MX|ST_1|Cluster_NA|2019/07/17`. The site comes from field 2 and the year from the text before the first `/` in field 4.

### Subsample targets

Hard-coded in `HIV_subsample_incidence.py`, derived from the lowest sequences-to-incidence ratio per year per location (DGE incidence; see Supplementary Table 1):

`ST_1` = CDMX, `ST_2` = MEX.

| Year | Location | Available | Incidence † | Seqs/incidence | Sampled |
|------|----------|----------:|------------:|---------------:|--------:|
| 2019 | CDMX | 571  | 14.2 | 40.21 | 310 |
| 2019 | MEX  | 229  | 10.5 | **21.81** ‡ | 229 |
| 2020 | CDMX | 705  | 7.6  | 92.76 | 166 |
| 2020 | MEX  | 254  | 5.1  | 49.80 | 111 |
| 2021 | CDMX | 1931 | 13.2 | 146.29 | 288 |
| 2021 | MEX  | 738  | 9.8  | 75.31 | 214 |
| 2022 | CDMX | 1899 | 14.4 | 131.88 | 314 |
| 2022 | MEX  | 751  | 12.1 | 62.07 | 264 |
| **Total** | | **7078** | | | **1896** |

† New HIV cases per 100,000 inhabitants per year, DGE (Greater Mexico City).
‡ Smallest sequences-to-incidence ratio; this cell sets the scale and is kept whole.

**Derivation.** Each target is `round(incidence × 21.81)`, where 21.81 is the smallest
sequences-to-incidence ratio across the eight cells (MEX 2019 = 229 / 10.5). All eight
hard-coded targets reproduce exactly under this rule, so the subsample is proportional
to reported incidence rather than to how much each location was sequenced. "Available"
is the full CIENI-CEC count for that year and location, and matches the supplementary table.

## Requirements

- Python 3
- Biopython (`pip install biopython`)

## Usage

```bash
python HIV_subsample_incidence.py   # input.fasta -> R1.fasta
python HIV_subsample_count.py       # prints the count table for R1.fasta
```

