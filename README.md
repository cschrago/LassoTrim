# LassoTrim_ScikitLearn

**LassoTrim_ScikitLearn** is a Python script designed to filter aligned nucleotide sequences by identifying phylogenetically informative sites using two strategies: LASSO-based site selection and Shannon entropy filtering. This tool is particularly suited for preprocessing alignments prior to phylogenomic inference.

## Features

- **LASSO Mode**: 
  - Performs replicate tree generation and site-wise log-likelihood calculations using **IQ-TREE 2**.
  - Applies **Lasso regression** via `scikit-learn` to select informative sites.
  - Outputs a filtered FASTA alignment including only the selected sites.

- **Entropy Mode**:
  - Calculates Shannon entropy for each alignment column.
  - Retains columns with entropy ≥ 0.5.
  - Produces a filtered FASTA file of high-entropy sites.

## Dependencies

- Python ≥ 3.6  
- Required Python packages:
  - `biopython`
  - `numpy`
  - `pandas`
  - `scikit-learn`
- External software:
  - `IQ-TREE 2` must be installed and available in the system path.

## Usage

```bash
python3 LassoTrim_ScikitLearn.py -s <alignment.fasta> [-l 1|0] [-r <replicates>]
```

### Arguments

| Argument        | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `-s`, `--seqfile`   | Input FASTA file with aligned sequences (required)                         |
| `-l`, `--lasso`     | Use LASSO analysis (`1`) or entropy filtering (`0`, default)              |
| `-r`, `--replicates`| Number of random tree replicates for LASSO site-wise likelihoods (default: 10000) |

## Output

Depending on the selected mode, the script generates the following:

### LASSO Mode (`-l 1`)
- `LassoClassification.txt`: Binary vector indicating informative sites (1 = selected).
- `<input>.lasso.info.fasta`: Filtered alignment with only LASSO-selected sites.
- `<input>.siteloglike.txt`: Site-wise log-likelihood matrix used in LASSO regression.

### Entropy Mode (`-l 0`)
- `<input>.entropy.info.fasta`: Filtered alignment including sites with entropy ≥ 0.5.

## Example

```bash
python3 LassoTrim_ScikitLearn.py -s example_aln.fasta -l 1 -r 5000
```

This command runs the LASSO pipeline on `example_aln.fasta`, generating 5000 replicate trees and producing a filtered alignment of informative sites.
