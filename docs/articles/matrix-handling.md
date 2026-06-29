# Matrix Handling

## Morphological matrix filtering

The function `filterMissing` deletes taxa and/or characters containing
only missing data. In the following example, the output file will be
saved as `test_filterMissing_FILTERED.nexus`:

``` r
library(RNODE)

filterMissing(input="../testdata/test_filterMissing.nexus", 
              input_format="nexus",
              output_path="../testdata/test_filterMissing",
              missing="both")
```

The function `filterInvariants` deletes invariant characters, which is
useful to accelerate graph searches. In Maximum Likelihood and Bayesian
analyses using the MKv model with ascertainment bias correction (ASC),
invariants must be deleted. Here, we follow the definition of invariant
from IQ-Tree, characterized by: (1) constant sites containing only a
single character state in all sequences, (2) partially constant sites (N
and/or -), and (3) ambiguously constant sites (e.g. C, Y and -). In the
following example, 122 invariants are detected.

``` r
filterInvariants(input="../testdata/015_MORPH_data.nexus",
                 input_format = "nexus",
                 output_index="../testdata/015_MORPH_data")
```

## Shared taxa filtering

The function `filterSharedTaxa` filters two molecular or morphological
matrices by taxon overlap. The parameter `level` controls whether the
output keeps only shared taxa, shared taxa plus taxa unique to input 1,
or shared taxa plus taxa unique to input 2. When a taxon is retained
from only one input, the missing row in the other matrix is filled with
`?`.

``` r
# Keep only taxa shared by the two matrices
filterSharedTaxa(input1="../testdata/013_MORPH_data.tnt", input1_format="tnt",
                 input2="../testdata/013_MOL_data.tnt", input2_format="tnt",
                 output_path="../testdata/013_TE_data2.tnt",
                 output_concatenate=TRUE,
                 level="strictly_shared")

# Keep shared taxa plus taxa unique to input1
filterSharedTaxa(input1="../testdata/013_MORPH_data.tnt", input1_format="tnt",
                 input2="../testdata/013_MOL_data.tnt", input2_format="tnt",
                 output_path="../testdata/013_TE_unique1.tnt",
                 output_concatenate=TRUE,
                 level="shared+unique1")

# Keep shared taxa plus taxa unique to input2
filterSharedTaxa(input1="../testdata/013_MORPH_data.tnt", input1_format="tnt",
                 input2="../testdata/013_MOL_data.tnt", input2_format="tnt",
                 output_path="../testdata/013_TE_unique2.tnt",
                 output_concatenate=TRUE,
                 level="shared+unique2")
```

If one of the input matrices contains protein data and the output is
TNT, `filterSharedTaxa` writes the TNT header `nstates 32`.

## Concatenation

The function `concatenate` merges two phylogenetic matrices column-wise.
Input can be NEXUS or TNT files (format is auto-detected from the file
extension) or matrices already loaded in R. Taxa present in only one of
the inputs are filled with `?` in the output. For instance, using the
data set from Whitcher et al. (2025), we concatenate a morphological
matrix (77 characters, 11 taxa) and a molecular matrix (45,165
characters, 28 taxa) into a total-evidence matrix (45,242 characters, 28
taxa):

``` r
concatenate(input1       = "../testdata/059_MORPH_data.tnt",
            input2       = "../testdata/059_MOL_data.tnt",
            output_file  = "../testdata/059_TE_data.tnt",
            output_format = "tnt")
```

The output file `059_TE_data.tnt` contains both character blocks
(`& [num]` for morphology, `& [num]` for molecular) separated by TNT
block headers. The 17 taxa absent from the morphological matrix are
automatically padded with `?` for the first 77 characters. The result
can also be captured as an R object without writing a file:

``` r
mat <- concatenate(input1 = "../testdata/059_MORPH_data.tnt",
                   input2 = "../testdata/059_MOL_data.tnt")
```

## Handling character states and orderings

The function `splitOrdFromUnord` splits a morphological matrix into
partitions of ordered and unordered characters based on a list of
ordered characters.

``` r
# Data input of list of ordered characters
list_ordered=c(1, 6, 7, 8, 10, 12, 13, 14, 17, 19, 23, 26, 31, 35, 41, 44, 45, 48, 51, 54, 55, 68, 71, 72, 92, 94, 96, 102, 105, 108, 109, 128, 129, 130, 131, 132, 135, 142, 144, 152, 153, 193)

splitOrdFromUnord(input="../testdata/048_MORPH_data.nex", 
                  input_format = "nexus",
                  output_index = "../testdata/048_MORPH", 
                  list_ordered=list_ordered)
```

The function `splitNoStates` splits characters from a morphological
matrix according to their number of character-states. This procedure has
been recommended to run phylogenetic analyses with the MK and MKv
models.

``` r
splitNoStates(input = "../testdata/015_MORPH_data.nexus", 
             input_format = "nexus", 
             output_index = "../testdata/015_MORPH_data", 
             ambiguity_addState = T, 
             inapplicable_addState = T, 
             log=T, 
             write=T)
```

## Format conversions

The functions `nexus2tnt` and `tnt2nexus` allow you to easily convert
matrices between NEXUS and TNT formats. They are specifically designed
to preserve data type headers (like `&[dna]` and `&[prot]`) and ordered
character indices (`ccode +`), which are often lost when using generic
tree-handling libraries. Ambiguities are also converted appropriately
(e.g. `(01)` to `[01]` in TNT).

``` r
# Convert NEXUS to TNT
nexus2tnt(input_file = "../testdata/048_MORPH_data.nex",
          output_file = "../testdata/048_converted.tnt")

# Convert TNT to NEXUS
tnt2nexus(input_file = "../testdata/059_MORPH_data.tnt",
          output_file = "../testdata/059_converted.nexus")
```
