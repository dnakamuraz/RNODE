# concatenate

Concatenate two or more phylogenetic matrices (morphological or
molecular) by columns. Inputs can be local NEXUS or TNT files, or
matrices already loaded in R. Taxa present in only one matrix are filled
with missing data ('?'). The input format and data types (morphology,
DNA, protein) are detected automatically if not explicitly defined in
headers. Also reports and translates ordered (additive) characters
safely to output matrix arrays.

## Usage

``` r
concatenate(
  input1 = NULL,
  input2 = NULL,
  input1_format = NULL,
  input2_format = NULL,
  output_file = NULL,
  output_format = NULL,
  multiple_input = NULL,
  gaps = NULL
)
```

## Arguments

- input1:

  First input: a local file path (NEXUS or TNT) or a matrix/data.frame
  already loaded in R. Ignored if \`multiple_input\` is provided.

- input2:

  Second input: a local file path (NEXUS or TNT) or a matrix/data.frame
  already loaded in R. Ignored if \`multiple_input\` is provided.

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

- multiple_input:

  A directory path containing multiple TNT/NEXUS files to concatenate.
  If provided, \`input1\` and \`input2\` are ignored.

- gaps:

  Logical controlling the gap-handling suffix written in TNT block
  headers (only relevant when \`output_format = 'tnt'\`). If \`TRUE\`,
  protein and DNA blocks are written as \`&\[prot gaps\]\` and \`&\[dna
  gaps\]\`. If \`FALSE\`, they are written as \`&\[prot nogaps\]\` and
  \`&\[dna nogaps\]\`. If \`NULL\` (default), no suffix is appended
  (\`&\[prot\]\` / \`&\[dna\]\`).

## Value

Invisibly returns the concatenated matrix.

## Author

Daniel YM Nakamura

## Examples

``` r
# Concatenate all files in a folder
concatenate(multiple_input = "testdata/partitioned/", output_file = "combined.tnt")
#> Error: Directory does not exist: testdata/partitioned/
```
