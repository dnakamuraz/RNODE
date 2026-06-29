# splitNoStates

Splits a morphological matrix according to the number of
character-states. This procedure has been recommended to run
phylogenetic analyses with the MK and MKv models (the 'K' refers to the
number of states). Khakurel et al. (2024) demonstrated that MK models
with high K values can underestimate the branch lengths, whereas MK
models with small K values can overstimate them. As such, some recent
studies have partitioned morphological characters according to their
number of states (e.g. Černý & Simonoff 2023).

## Usage

``` r
splitNoStates(
  input,
  input_format = NULL,
  ambiguity_addState = FALSE,
  inapplicable_addState = FALSE,
  output_index = "output",
  log = FALSE,
  write = TRUE
)
```

## Arguments

- input:

  Input file (morphological matrix, either loaded previously in R or
  loaded from a local file).

- ambiguity_addState:

  If T, ambiguities are counted as additional character-states (default:
  ambiguity_addState = F).

- inapplicable_addState:

  If T, inapplicable states are counted toward the sum of unique
  character-states.

- output_index:

  Output index (e.g. if the user specify it as "Desktop/Index", the
  output files will be "Desktop/Index_ORDERED.nexus" and
  "Desktop/Index_UNORDERED.nexus")

- log:

  If T, write a file locally reporting the destination of each
  character.

- write:

  If T, write the files locally.

- inpu_format:

  To load from a local file: 'nexus' or 'tnt'.

## References

Černý, D., & Simonoff, A. L. (2023). Statistical evaluation of character
support reveals the instability of higher-level dinosaur phylogeny.
Scientific Reports, 13(1), 9273.

Khakurel, B., Grigsby, C., Tran, T. D., Zariwala, J., Höhna, S., &
Wright, A. M. (2024). The fundamental role of character coding in
Bayesian morphological phylogenetics. Systematic biology, 73(5),
861-871.

## Author

Daniel YM Nakamura

## Examples

``` r
# Synthetic example
splitNoStates(input = "../testdata/015_MORPH_data.nexus", input_format = "nexus", output_index = "../testdata/015_MORPH_data", ambiguity_addState = T, inapplicable_addState = T, log=T, write=T)
#> Error in .UTFLines(filepath, encoding): File '../testdata/015_MORPH_data.nexus' not found.
```
