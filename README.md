# RNODE: Facilitating pre- and postprocessing of phylogenetic analyses

[![language](https://img.shields.io/badge/language-R-blue?style=flat&logo=r&logoColor=white)](https://www.r-project.org)
[![author](https://img.shields.io/badge/author-DYM_Nakamura-blue?logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=c0W8Cm8AAAAJ&hl=en)
[![license](https://img.shields.io/badge/license-GPL_v3-blue?logo=gnu&logoColor=white)](https://www.gnu.org/licenses/gpl-3.0.html)

**RNODE** is an R package to facilitate pre- and postprocessing of phylogenetic analyses, including (1) comparisons of topologies, branch lengths, support values, (2) comparison and consensus of DNA sequences, (3) manipulation of cladistic matrices, and (4) manipulation of trees.

Copyright (C) Daniel Y. M. Nakamura 2026

## Cite

If you use **RNODE**, please cite this repository.

## Installation

**RNODE** was tested in R. v. 4.5.2 and can be installed with the following command:

```
devtools::install_github("dnakamuraz/RNODE")
```

## Usage

For detailed tutorials and examples on how to use these functions, please refer to the package vignettes available on the website's Articles section:

- [Tree Comparisons](https://dnakamuraz.github.io/RNODE/articles/tree-comparisons.html): Identifying shared and unique clades, support comparisons, logistic regressions, branch lengths, and topological distances.
- [Matrix Handling](https://dnakamuraz.github.io/RNODE/articles/matrix-handling.html): Morphological filtering, shared taxa filtering, concatenating matrices, and formatting conversions (TNT / NEXUS).
- [Tree Manipulation](https://dnakamuraz.github.io/RNODE/articles/tree-manipulation.html): Mapping support values and branch lengths across different trees.

### Consensus sequences

`consensusFAS()` builds a consensus sequence from a FASTA alignment using IUPAC ambiguity codes and appends it to the alignment. Three `consensus` methods are available:

| `consensus =` | Rule | Example output |
|---|---|---|
| `"relaxed"` (default) | IUPAC union of all observed nucleotides; gaps/`N` ignored when real bases present | `ATCGTTAA` |
| `"conservative"` | Same, but positions where >50 % of sequences are gaps/`N` become `?` | `ATCGTT??` |
| `"both"` | Appends both sequences (`Consensus_relaxed` + `Consensus_conservative`) | — |

**`N` / gap handling** per column:

| Column composition | Relaxed | Conservative |
|---|---|---|
| All `N` | `N` | `N` |
| Mix of `N` and `-` | `N` | `?` |
| All `-` / `?` | `-` | `-` |
| Real nucleotide(s) present | IUPAC symbol | IUPAC symbol (or `?` if majority missing) |

```r
# Relaxed consensus (default)
result <- consensusFAS(input = "alignment.fasta")

# Conservative consensus, written to file
result <- consensusFAS(input = "alignment.fasta",
                       consensus = "conservative",
                       output = "consensus_alignment.fasta")

# Both methods appended in one call
result <- consensusFAS(input = "alignment.fasta",
                       consensus = "both",
                       output = "consensus_alignment.fasta")

# Align unaligned sequences first (requires MUSCLE in PATH)
result <- consensusFAS(input = "unaligned.fasta", align = TRUE)
```

