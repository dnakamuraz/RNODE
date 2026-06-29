# retrodictNodes

Given two trees with support values, creates a dataframe containing
support values of tree 1 and the occurrence of the clade in tree 2. This
dataframe can be used for logistic regressions testing whether support
values of tree 1 retrodict the occurrence of clades in tree 2.

## Usage

``` r
retrodictNodes(tree1, tree2, outgroup = NULL, root = NULL, dataframe = F)
```

## Arguments

- tree1:

  A .phylo tree that can be loaded using ape::read.tree for NEWICK files
  or TreeTools::ReadTntTree for TNT files. The .phylo must contain
  \$node.label.

- tree2:

  Another .phylo tree.

- outgroup:

  Optional. Specify outgroup taxa to remove (by default, the function
  assumes that the user does not want to remove outgroup taxa)

- root:

  Optional. Specify the same root for both trees, which is recommended
  to facilitate tree comparisons (by default, the function assumes that
  trees share the same root)

- dataframe:

  Optional. Write a TSV file in current directory containing the output
  dataframe (by default, no .TSV is written).

## Author

Daniel YM Nakamura

## Examples

``` r
# Example 1 (simplest case)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
retrodictNodes (tree1, tree2)
#> All required parameters provided.
#> All required parameters provided.
#> [1] ""
#> [1] "Tree comparisons done!"
#> [1] "Number of shared clades:  3"
#> [1] "Tree 1: Total number of clades = 3 ; Mean support = 53.5"
#> [1] "Support of shared clades in tree 1: 32\342\200\22375 (53.5)"
#> [1] "Tree 2: Total number of clades = 3 ; Mean support = 50"
#> [1] "Support of shared clades in tree 2: 47\342\200\22353 (50)"
#> All required parameters provided.
#> [1] "Both trees with support values."
#> character(0)
#> character(0)
#>   support_tree1 occurrence_tree2
#> 2            32                1
#> 3            75                1
```
