# uniqueNodes

Creates two separate dataframes containing the list of descendants and
support values of unique clades between two trees.

## Usage

``` r
uniqueNodes(
  tree1,
  tree2,
  composition = T,
  outgroup = NULL,
  root = NULL,
  dataframe = F,
  dataframe1.name = "Tree1_unique.clades.tsv",
  dataframe2.name = "Tree2_unique.clades.tsv",
  plotTrees = F,
  node.numbers = T,
  tree.width = 10,
  tree.height = 10,
  tree.fsize = 0.5,
  tree.adj = c(-1.5, 0.5),
  tree.cex = 2,
  sup = T,
  sup.adj1 = c(-1, 1),
  sup.adj2 = c(1, 1),
  sup.cex = 1,
  output.tree = "trees_unique_nodes.pdf"
)
```

## Arguments

- tree1:

  A .phylo tree that can be loaded using ape::read.tree for NEWICK files
  or TreeTools::ReadTntTree for TNT files

- tree2:

  Another .phylo tree

- composition:

  Optional. Specify if composition of corresponding clades should be
  present in the dataframe (by default, composition = T)

- outgroup:

  Optional. Specify outgroup taxa to remove (by default, outgroup = F
  assumes that the user does not want to remove outgroup taxa)

- root:

  Optional. Specify the same root for both trees, which is recommended
  to facilitate tree comparisons (by default, root = F assumes that
  trees share the same root)

- dataframe:

  Optional. Write a TSV file in current directory containing the output
  dataframe (by default, dataset = T).

- dataframe1.name:

  Optional. Name to write the dataframe 1.

- dataframe2.name:

  Optional. Name to write the dataframe 2.

- plotTrees:

  Optional. Plot the two trees after taxa pruning in `PDF` format. If
  `plot = T`, the user should also adjust `PDF` dimensions (e.g.
  `width = 8`, `height = 8`), label size (e.g. `fsize = 4`), and
  position and size of support values (e.g. `adj = c(-1.5,0.5)`,
  `cex = 0.6`).

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

- sup:

  Optional. Optional. Plot support values (by default, sup = T).

- sup.adj1:

  Optional. If sup = T, adjust horizontal and vertical position of
  support values on tree1.

- sup.adj2:

  Optional. If sup = T, adjust horizontal and vertical position of
  support values on tree2.

- sup.cex:

  Optional. If sup = T, adjust font size of support values on tree1.

- output.tree:

  Optional. The output file name of tree 1 if plotTrees = T.

## Author

Daniel YM Nakamura, Taran Grant

## Examples

``` r
# Example 1 (identify unique nodes)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
uniqueNodes (tree1, tree2)
#> All required parameters provided.
#> [1] "Both trees with support values."
#> character(0)
#> character(0)
#> [[1]]
#> [1] Node    Support
#> <0 rows> (or 0-length row.names)
#> 
#> [[2]]
#> [1] Node    Support
#> <0 rows> (or 0-length row.names)
#> 

# Example 2 (count unique nodes)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
nrow (uniqueNodes (tree1, tree2, composition = F))
#> All required parameters provided.
#> [1] "Both trees with support values."
#> character(0)
#> character(0)
#> NULL

# Exmple 3 (highlight unique nodes in plotting)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
uniqueNodes (tree1, tree2, plotTrees=T, node.numbers=F, tree.width=14, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2)
#> All required parameters provided.
#> [1] "Both trees with support values."
#> character(0)
#> character(0)
#> [[1]]
#> [1] Node    Support
#> <0 rows> (or 0-length row.names)
#> 
#> [[2]]
#> [1] Node    Support
#> <0 rows> (or 0-length row.names)
#> 
```
