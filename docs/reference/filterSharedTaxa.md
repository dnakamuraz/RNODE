# filterSharedTaxa

Filter two molecular or morphological matrices by shared and unique taxa
(rows). Can output concatenated or separated matrices in NEXUS or TNT
format.

## Usage

``` r
filterSharedTaxa(
  input1,
  input1_format = NULL,
  input2,
  input2_format = NULL,
  output_path = NULL,
  output_format = NULL,
  output_concatenate = TRUE,
  return_as_list = FALSE,
  level = "strictly_shared"
)
```

## Arguments

- input1:

  First input file (molecular or morphological matrix in nexus, fasta,
  or tnt format) or already loaded as a matrix object in R.

- input1_format:

  To load from a local file: 'nexus', 'fasta', or 'tnt'. Defaults to
  NULL (auto-detects).

- input2:

  Second input file (molecular or morphological matrix in nexus, fasta,
  or tnt format) or already loaded as a matrix object in R.

- input2_format:

  To load from a local file: 'nexus', 'fasta', or 'tnt'. Defaults to
  NULL (auto-detects).

- output_path:

  Output path. For NEXUS output, "Desktop/Index" will output
  "Desktop/Index_SHARED.nexus" (concatenated) or
  "Desktop/Index_SHARED_1.nexus" and "Desktop/Index_SHARED_2.nexus"
  (separated). For TNT output, "Desktop/Index" will output
  "Desktop/Index.tnt" (concatenated) or "Desktop/Index_SHARED_1.tnt" and
  "Desktop/Index_SHARED_2.tnt" (separated). If NULL, returns matrices as
  a list without writing files.

- output_format:

  Output file format: 'nexus' or 'tnt'. If NULL (default), the format is
  inferred from output_path: paths ending in .tnt write TNT; all others
  write NEXUS.

- output_concatenate:

  Logical. If TRUE (default), concatenates the filtered matrices. If
  FALSE, returns them separately.

- return_as_list:

  Logical. If TRUE, returns the result as a list of matrices instead of
  writing to file. Default is FALSE.

- level:

  Taxon filtering level. Use 'strictly_shared' to keep only taxa shared
  by input1 and input2, 'shared+unique1' to keep shared taxa plus taxa
  unique to input1, or 'shared+unique2' to keep shared taxa plus taxa
  unique to input2. Missing data for taxa absent from one matrix are
  filled with '?'.

## Author

Daniel YM Nakamura

## Examples

``` r
# Example with concatenated NEXUS output
filterSharedTaxa(
  input1 = "testdata/file1.nexus", input1_format = "nexus",
  input2 = "testdata/file2.nexus", input2_format = "nexus",
  output_path = "testdata/shared", output_concatenate = TRUE
)
#> Warning: cannot open file 'testdata/file1.nexus': No such file or directory
#> Error in file(con, "r"): cannot open the connection

# Example with concatenated TNT output (protein data)
filterSharedTaxa(
  input1 = "testdata/013_MORPH_data.tnt", input1_format = "tnt",
  input2 = "testdata/013_MOL_data.tnt", input2_format = "tnt",
  output_path = "testdata/013_TE_data", output_format = "tnt",
  output_concatenate = TRUE, level = "strictly_shared"
)
#> Warning: cannot open file 'testdata/013_MORPH_data.tnt': No such file or directory
#> Error in file(con, "r"): cannot open the connection

# Example with separated output
filterSharedTaxa(
  input1 = "testdata/file1.nexus", input1_format = "nexus",
  input2 = "testdata/file2.nexus", input2_format = "nexus",
  output_path = "testdata/shared", output_concatenate = FALSE
)
#> Warning: cannot open file 'testdata/file1.nexus': No such file or directory
#> Error in file(con, "r"): cannot open the connection
```
