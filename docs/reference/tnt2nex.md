# tnt2nex

Converts a TNT matrix into NEXUS format while preserving header data
type and ordered character structures.

## Usage

``` r
tnt2nex(input_file, output_file, iupac_ambiguity = FALSE)
```

## Arguments

- input_file:

  Path to the input TNT file.

- output_file:

  Path to the output NEXUS file.

- iupac_ambiguity:

  Logical; if TRUE, uses IUPAC symbols for ambiguities (e.g. R, Y),
  otherwise outputs between parentheses (e.g. (CT)).

## Author

Daniel YM Nakamura
