# filterTaxaList

Filter a molecular or morphological matrix by keeping or deleting
specific taxa based on a provided list.

## Usage

``` r
filterTaxaList(
  input,
  input_format = NULL,
  list,
  filter = "keep",
  output = NULL,
  output_format = NULL
)
```

## Arguments

- input:

  A local file path (molecular or morphological matrix in nexus or tnt
  format) or already loaded as a matrix object in R.

- input_format:

  To load from a local file: 'nexus', 'tnt', or NULL. Defaults to NULL
  (auto-detects from extension).

- list:

  A local .txt file path or an R object (vector, list, or data.frame)
  containing the list of terminals, with one terminal name per
  row/element.

- filter:

  Character string. Use 'delete' to remove taxa contained in the list,
  or 'keep' to retain only the taxa contained in the list.

- output:

  Output path (e.g., "Desktop/Filtered"). If NULL, returns the matrix as
  an R object without writing files.

- output_format:

  Output file format: 'nexus', 'tnt', or NULL. If NULL, the format is
  inferred from the \`output\` extension, or simply returns as an R
  object if \`output\` is also NULL.

## Author

Daniel YM Nakamura

## Examples

``` r
if (FALSE) { # \dontrun{
# Keep only the taxa specified in a text file
filterTaxaList(
    input = "testdata/matrix.tnt",
    list = "testdata/taxa_to_keep.txt",
    filter = "keep",
    output = "testdata/filtered_matrix.tnt"
)

# Delete taxa using an R vector and return as an object
taxa_to_remove <- c("TaxonA", "TaxonB")
filtered_mat <- filterTaxaList(
    input = "testdata/matrix.nexus",
    list = taxa_to_remove,
    filter = "delete"
)
} # }
```
