# splitOrdFromUnord

Splits a morphological matrix containing ordered and unordered
characters into two matrices (ordered and unordered). This is useful to
run phylogenetic analyses with MK and ORDERED models in IQTREE,
especially if a concatenated matrix containing invariant characters is
given as input.

## Usage

``` r
splitOrdFromUnord(input, input_format, output_index, list_ordered)
```

## Arguments

- input:

  Input file (concatenated morphological matrix loaded locally or
  already loaded in R).

- input_format:

  To load from a local file: 'nexus' or 'tnt'.

- output_index:

  Output index (e.g. if the user specify it as "Desktop/Index", the
  output files will be "Desktop/Index_ORDERED.nexus" and
  "Desktop/Index_UNORDERED.nexus")

- list_ordered:

  List of ordered characters e.g. c(1, 3, 9, 13). Character numbering
  starts with 1 (even if input data is .TNT).

## Author

Daniel YM Nakamura

## Examples

``` r
# Example
splitOrdFromUnord(input="../testdata/048_MORPH_data.nex", output_index = "../testdata/048_MORPH", list_ordered=list_ordered)
#> Error in splitOrdFromUnord(input = "../testdata/048_MORPH_data.nex", output_index = "../testdata/048_MORPH",     list_ordered = list_ordered): argument "input_format" is missing, with no default
```
