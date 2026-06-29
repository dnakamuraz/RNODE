# compareBranchLength

Given two trees with branch lengths, return a dataframe containing
lengths from matching branches.

## Usage

``` r
compareBranchLength(
  tree1,
  tree2,
  write = NULL,
  outgroup = NULL,
  root = NULL,
  composition = F,
  unique = F
)
```

## Arguments

- tree1:

  A phylo object without branch lengths

- tree2:

  A phylo object with branch lengths

- write:

  Optional. Specify the name of the dataframe file to be written locally
  (nothing is written if this parameter is not specified)

- outgroup:

  Optional. Specify outgroup taxa to remove (by default, outgroup = F
  assumes that the user does not want to remove outgroup taxa)

- root:

  Optional. Specify the same root for both trees, which is recommended
  to facilitate tree comparisons (by default, root = F assumes that
  trees share the same root)

- composition:

  Optional. Specify if composition of corresponding clades should be
  present in the dataframe (by default, composition = F)

- unique:

  Optional. Show unique clades (default: unique = F)

## Author

Daniel YM Nakamura

## Examples

``` r
# Example 1 (identify unique nodes)
tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
compareBranchLength (tree1, tree2)
#> All required parameters provided.
#> Error in data.frame(ParentNode = tree$edge[, 1], ChildNode = tree$edge[,     2], EdgeType = edge_type, Clade = sapply(clades, paste, collapse = ","),     EdgeLength = tree$edge.length, stringsAsFactors = FALSE): arguments imply differing number of rows: 8, 0
```
