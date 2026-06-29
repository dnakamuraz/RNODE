# renameTaxa

Renames terminals (taxa) in a phylogenetic matrix (TNT, NEXUS, or an R
matrix) using a provided synonym dictionary. It can optionally delete
taxa that are not specified in the synonym file and output the resulting
matrix to a local file.

## Usage

``` r
renameTaxa(
  input,
  input_format = NULL,
  synonym_file,
  delete_missing = FALSE,
  output = NULL,
  output_format = NULL
)
```

## Arguments

- input:

  A local file path (NEXUS or TNT) or a matrix/data.frame already loaded
  in R.

- input_format:

  Format of the input file: 'nexus' or 'tnt'. Detected automatically if
  NULL (default).

- synonym_file:

  Path to a text file containing the synonym list. In each row, the
  first string is the new name, and subsequent strings (separated by
  spaces or tabs) are the old names/synonyms.

- delete_missing:

  Logical. If TRUE, deletes all terminals in the input matrix that are
  not present in the \`synonym_file\`. Default is FALSE.

- output:

  Output file path (e.g. "output.tnt"). If NULL, the result is returned
  as an R object only.

- output_format:

  Output format: 'nexus' or 'tnt'. If NULL, inferred from the extension
  of \`output\`. Ignored when output is NULL.

## Value

Invisibly returns the matrix with renamed (and optionally filtered)
terminals.

## Author

\[Your Name\]

## Examples

``` r
if (FALSE) { # \dontrun{
# Rename taxa in a TNT file and write the result back to TNT
renameTaxa(
    input = "input.tnt",
    synonym_file = "synonyms.txt",
    delete_missing = TRUE,
    output = "renamed_output.tnt"
)

# Read from a loaded R matrix and return an R object only
mat <- renameTaxa(
    input = my_r_matrix,
    synonym_file = "synonyms.txt"
)
} # }
```
