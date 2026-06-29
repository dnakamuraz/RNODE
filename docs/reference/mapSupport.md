# mapSupport

Given a tree 1 without support values (e.g. strict consensus of optimal
trees) and tree 2 with support values (e.g. majority consensus from
suboptimal bootstrap pseudo-replicates), return tree 1 with support
values from shared clades with tree 2.

## Usage

``` r
mapSupport(
  tree1,
  tree2,
  write = NULL,
  outgroup = NULL,
  root = NULL,
  plotTrees = F,
  node.numbers = T,
  tree.width = 10,
  tree.height = 10,
  tree.fsize = 0.5,
  tree.adj = c(-1.5, 0.5),
  tree.cex = 2,
  tree.output = "trees_unique_nodes.pdf"
)
```

## Arguments

- tree1:

  A phylo object without support values

- tree2:

  A phylo object with support values

- write:

  Optional. Specify the name of tree file to be written locally (nothing
  is written if this parameter is not specified)

- outgroup:

  Optional. Specify outgroup taxa to remove (by default, outgroup = F
  assumes that the user does not want to remove outgroup taxa)

- root:

  Optional. Specify the same root for both trees, which is recommended
  to facilitate tree comparisons (by default, root = F assumes that
  trees share the same root)

- plotTrees:

  Optional. Plot the two trees (tree1 with support values on left and
  tree2 on right) in `PDF` format. If `plot = T`, the user should also
  adjust `PDF` dimensions (e.g. `width = 8`, `height = 8`), label size
  (e.g. `fsize = 4`), and position and size of support values (e.g.
  `adj = c(-1.5,0.5)`, `cex = 0.6`).

- node.numbers:

  Optional. If plotTrees = T, show node index (do not confuse with
  support values'by default, True).

- tree.width:

  Optional. Width of trees in PDF if plotTrees = T.

- tree.height:

  Optional. Height of trees in PDF if plotTrees = T.

- tree.fsize:

  Optional. Font size in PDF if plotTrees = T.

- tree.adj:

  Optional. Adjust horizontal and vertical position if plotTrees = T.

- tree.cex:

  Optional. Adjust support size in nodes if plotTrees = T.

- tree.output:

  Optional. Name of the output figure.

## Author

Daniel YM Nakamura

## Examples

``` r
# Example 1 (identify unique nodes)
tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
mapSupport (tree1, tree2)
#> All required parameters provided.
#> [[1]]
#> 
#> Phylogenetic tree with 5 tips and 4 internal nodes.
#> 
#> Tip labels:
#>   t1, t3, t2, t4, t5
#> Node labels:
#>   NA, 94, NA, 47
#> 
#> Rooted; no branch length.
#> 
#> [[2]]
#>      tr1 tr2 tr2_support
#> [1,]   6   6          NA
#> [2,]   7   7          94
#> [3,]   8  NA          NA
#> [4,]   9   9          47
#> 
```
