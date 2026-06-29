# concatenate

Concatenate two phylogenetic matrices (morphological or molecular) by
columns. Inputs can be local NEXUS or TNT files, or matrices already
loaded in R. Taxa present in only one matrix are filled with missing
data ('?'). The input format and data types (morphology, DNA, protein)
are detected automatically if not explicitly defined in headers. Also
reports and translates ordered (additive) characters safely to output
matrix arrays.

## Usage

``` r
concatenate(
  input1,
  input2,
  input1_format = NULL,
  input2_format = NULL,
  output_file = NULL,
  output_format = NULL
)
```

## Arguments

- input1:

  First input: a local file path (NEXUS or TNT) or a matrix/data.frame
  already loaded in R.

- input2:

  Second input: a local file path (NEXUS or TNT) or a matrix/data.frame
  already loaded in R.

- input1_format:

  Format of input1 file: 'nexus' or 'tnt'. Detected automatically if
  NULL (default).

- input2_format:

  Format of input2 file: 'nexus' or 'tnt'. Detected automatically if
  NULL (default).

- output_file:

  Output file path (e.g. "testdata/059_TE_data.tnt"). If NULL, the
  result is returned as an R object only.

- output_format:

  Output format: 'nexus' or 'tnt'. If NULL, inferred from the extension
  of output_file ('.tnt' -\> tnt, otherwise nexus). Ignored when
  output_file is NULL.

## Value

Invisibly returns the concatenated matrix.

## Author

Daniel YM Nakamura

## Examples

``` r
# Concatenate two TNT files
concatenate(
  input1 = "testdata/059_MORPH_data.tnt",
  input2 = "testdata/059_MOL_data.tnt",
  output_file = "testdata/059_TE_data.tnt",
  output_format = "tnt"
)
#> Warning: cannot open file 'testdata/059_MORPH_data.tnt': No such file or directory
#> Error in file(con, "r"): cannot open the connection

# Return as R object without writing
mat <- concatenate(
  input1 = "testdata/059_MORPH_data.tnt",
  input2 = "testdata/059_MOL_data.tnt"
)
#> Warning: cannot open file 'testdata/059_MORPH_data.tnt': No such file or directory
#> Error in file(con, "r"): cannot open the connection
```
