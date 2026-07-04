# tnt2fas

Converts a TNT matrix into FASTA format. Only works for DNA and PROTEIN
data.

## Usage

``` r
tnt2fas(input_file, output_file, iupac_ambiguity = TRUE)
```

## Arguments

- input_file:

  Path to the input TNT file.

- output_file:

  Path to the output FASTA file.

- iupac_ambiguity:

  Logical; if TRUE, uses IUPAC symbols for ambiguities (default TRUE for
  FASTA).

## Author

Daniel YM Nakamura
