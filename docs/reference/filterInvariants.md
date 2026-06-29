# filterInvariants

Delete invariant characters in a morphological matrix. We follow the
definition of invariant site from IQ-Tree: (1) constant sites containing
only a single character state in all sequences, (2) partially constant
sites (N and/or -), and (3) ambiguously constant sites (e.g. C, Y and
-).

## Usage

``` r
filterInvariants(
  input,
  input_format = NULL,
  output_format = NULL,
  output_index = NULL
)
```

## Arguments

- input:

  Input file (molecular or morphological matrix loaded locally or
  already loaded as a matrix object in R).

- input_format:

  To load from a local file: 'nexus', 'tnt', or NULL for auto-detect.
  Default is NULL.

- output_format:

  To write a local file: 'tnt', 'nexus', or NULL (default; only returns
  as an R object).

- output_index:

  Output index (e.g. if the user specify it as "Desktop/Index", the
  output files will be "Desktop/Index_onlyVARIANTS.nexus"). Required if
  output_format is not NULL.

## Author

Daniel YM Nakamura

## Examples

``` r
# Example returning only R object:
mat_filtered <- filterInvariants(input = "../testdata/015_MORPH_data.nexus")
#> Warning: cannot open file '../testdata/015_MORPH_data.nexus': No such file or directory
#> Error in file(con, "r"): cannot open the connection

# Example writing locally as TNT:
filterInvariants(input = "../testdata/015_MORPH_data.nexus", output_format = "tnt", output_index = "../testdata/015_MORPH_data")
#> Warning: cannot open file '../testdata/015_MORPH_data.nexus': No such file or directory
#> Error in file(con, "r"): cannot open the connection
```
