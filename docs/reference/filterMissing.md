# filterMissing

Filter taxa (rows) and/or characters (columns) in a matrix containing
only missing data.

## Usage

``` r
filterMissing(
  input,
  input_format = NULL,
  output,
  output_format = "nexus",
  missing = "row"
)
```

## Arguments

- input:

  Input file (morphological matrix in Nexus format).

- input_format:

  To load from a local file: 'nexus' or 'tnt'.

- output:

  Output path (e.g. if the user specify it as "Desktop/Index", the
  output files will be "Desktop/Index_ORDERED.nexus" and
  "Desktop/Index_UNORDERED.nexus")

- output_format:

  Format to write the output matrix. Options: "nexus" (default) or
  "tnt".

- missing:

  Parameter specifying if rows and/or columns in which all cells are
  missing data (?) should be removed. Options: "row" (default i.e.
  terminals), "column" (i.e. characters, transformation series), "both".

## Author

Daniel YM Nakamura

## Examples

``` r
if (FALSE) { # \dontrun{
# Example
filterMissing(input = "testdata/test_filterMissing.nexus", output = "testdata/test", missing = "row")
} # }
```
