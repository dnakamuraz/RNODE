# sharedNodes

Compare support values of shared clades between two trees. The outputs
are (1) basic statistics about number of shared clades and support
values; (2) a dataframe with node labels, descendants, and support
values of shared clades, which facilitates descriptive and statistical
comparisons of clade composition and support between corresponding
nodes.

## Usage

``` r
sharedNodes(
  tree1,
  tree2,
  composition = F,
  outgroup = NULL,
  root = NULL,
  plotTrees = F,
  node.numbers = T,
  tree.width = NULL,
  tree.height = NULL,
  tree.fsize = NULL,
  tree.adj = NULL,
  tree.cex = NULL,
  tree2.direction = "rightwards",
  output.tree1 = "tree1_pruned.pdf",
  output.tree2 = "tree2_pruned.pdf",
  tanglegram = F,
  output.tangletree = "tanglegram_comparison.pdf",
  tanglegram.margin = 7,
  tanglegram.lab.cex = 0.5,
  tanglegram.edge.lwd = c(0.1, 0.1),
  tanglegram.lwd = 1,
  tanglegram.axes = F,
  tanglegram.width = 5,
  tanglegram.height = 5,
  tanglegram.colors = T,
  dataframe = F,
  output.dataframe = "shared.nodes.tsv",
  messages = T,
  spearman = F,
  write.pruned = F,
  write.pruned1.name = "pruned_tree1.nwk",
  write.pruned2.name = "pruned_tree2.nwk"
)
```

## Arguments

- tree1:

  A `phylo` tree that can be loaded using
  [`ape::read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html) for
  `NEWICK` files or
  [`TreeTools::ReadTntTree`](https://ms609.github.io/TreeTools/reference/ReadTntTree.html)
  for `TNT` files. The `phylo` must contain the element `$node.label`.

- tree2:

  Another `phylo` tree

- composition:

  Optional. Specify if composition of corresponding clades should be
  present in the dataframe (by default, `composition = F`)

- outgroup:

  Optional. Specify outgroup taxa to remove (by default, the function
  assumes that the user does not want to remove outgroup taxa; i.e.
  `outgroup = NULL`)

- root:

  Optional. Specify the same root for both trees, which is recommended
  to facilitate tree comparisons (by default, the function assumes that
  trees share the same root; i.e. `root = NULL`)

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

- tree2.direction:

  Optional. Adjust the direction of tree2 ("leftwards" vs "rightwards")
  if plotTrees = T.

- output.tree1:

  Optional. The output file name of tree 1 if plotTrees = T.

- output.tree2:

  Optional. The output file name of tree 2 if plotTrees = T.

- tanglegram:

  Optional. Plot a tanglegram minimizing the number of crosses of lines
  linking two trees in `PDF` format. If the input tree has no branch
  length, uniform lengths are simulated to enable visualization.

- output.tangletree:

  Optional. The output file name if tanglegram = T.

- tanglegram.margin:

  Optional. Distance between tangle trees (default = 7) if tanglegram =
  T.

- tanglegram.lab.cex:

  Optional. Size of leaf names (default = 0.5) if tanglegram = T.

- tanglegram.edge.lwd:

  Optional. Thickness of edges of tangle trees (default = c(0.1, 0.1))
  if tanglegram = T.

- tanglegram.lwd:

  Optional. Thickness of edges connecting both tangle trees (default
  = 1) if tanglegram = T.

- tanglegram.axes:

  Optional. Show scale in tangle trees (default = F) if tanglegram = T.

- tanglegram.width:

  Optional. Width of tangle trees in PDF (default = 5) if tanglegram =
  T.

- tanglegram.height:

  Optional. Height of tangle trees in PDF (default = 5) if tanglegram =
  T.

- tanglegram.colors:

  Optional. Show color in edges connecting both tangle trees (default
  = T) if tanglegram = T.

- dataframe:

  Optional. Write a `TSV` file in current directory containing the
  output dataframe (by default, no `TSV` is written).

- output.dataframe:

  Optional. The output file name of the dataframe if dataframe = T.

- spearman:

  Optional. Test the correlation between support values using a Spearman
  test (by default, `spearman = F`).

- write.pruned:

  Write the pruned trees.

- write.pruned1.name:

  Output file name of pruned tree1 if write.pruned = T.

- write.pruned2.name:

  Output file name of pruned tree2 if write.pruned = T.

## Author

Daniel YM Nakamura, Taran Grant

## Examples

``` r
# Example 1 (simplest case)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
sharedNodes (tree1, tree2)
#> All required parameters provided.
#> [1] ""
#> [1] "Tree comparisons done!"
#> [1] "Number of shared clades:  3"
#> [1] "Tree 1: Total number of clades = 3 ; Mean support = 53.5"
#> [1] "Support of shared clades in tree 1: 32\342\200\22375 (53.5)"
#> [1] "Tree 2: Total number of clades = 3 ; Mean support = 50"
#> [1] "Support of shared clades in tree 2: 47\342\200\22353 (50)"
#>   Node_Tree_1 Node_Tree_2 Support_Tree_1 Support_Tree_2
#> 1           5           5                              
#> 2           6           6             32             53
#> 3           7           7             75             47

# Example 2 (show internal topology of each node, remove outgroup taxa t9 and t8, and reroot in t1 in both trees)
sharedNodes (tree1, tree2, composition=T, outgroup=c("t9", "t8"), root="t1",)
#> All required parameters provided.
#> [1] ""
#> [1] "Tree comparisons done!"
#> [1] "Number of shared clades:  3"
#> [1] "Tree 1: Total number of clades = 3 ; Mean support = 53.5"
#> [1] "Support of shared clades in tree 1: 32\342\200\22375 (53.5)"
#> [1] "Tree 2: Total number of clades = 3 ; Mean support = 50"
#> [1] "Support of shared clades in tree 2: 47\342\200\22353 (50)"
#>   Node_Tree_1 Node_Tree_2 Support_Tree_1 Support_Tree_2 Composition_Tree_1
#> 1           5           5                               (((t4,t5),t3),t1);
#> 2           6           6             32             53      ((t4,t5),t3);
#> 3           7           7             75             47           (t4,t5);
#>   Composition_Tree_2
#> 1 (((t4,t5),t3),t1);
#> 2      ((t4,t5),t3);
#> 3           (t4,t5);

# Example 3 (plot  two trees)
sharedNodes (tree1, tree2, plotTrees=T, tree.width=8, tree.height=8, tree.fsize=3, tree.adj=c(-1.5,0.5), tree.cex=0.6)
#> All required parameters provided.
#> [1] ""
#> [1] "Tree comparisons done!"
#> [1] "Number of shared clades:  3"
#> [1] "Tree 1: Total number of clades = 3 ; Mean support = 53.5"
#> [1] "Support of shared clades in tree 1: 32\342\200\22375 (53.5)"
#> [1] "Tree 2: Total number of clades = 3 ; Mean support = 50"
#> [1] "Support of shared clades in tree 2: 47\342\200\22353 (50)"
#>   Node_Tree_1 Node_Tree_2 Support_Tree_1 Support_Tree_2
#> 1           5           5                              
#> 2           6           6             32             53
#> 3           7           7             75             47

# Example 4 (plot tangle trees)
sharedNodes (tree1, tree2, tanglegram=T, tanglegram.margin=7, tanglegram.lab.cex=0.5)
#> All required parameters provided.
#> [1] ""
#> [1] "Tree comparisons done!"
#> [1] "Number of shared clades:  3"
#> [1] "Tree 1: Total number of clades = 3 ; Mean support = 53.5"
#> [1] "Support of shared clades in tree 1: 32\342\200\22375 (53.5)"
#> [1] "Tree 2: Total number of clades = 3 ; Mean support = 50"
#> [1] "Support of shared clades in tree 2: 47\342\200\22353 (50)"
#>   Node_Tree_1 Node_Tree_2 Support_Tree_1 Support_Tree_2
#> 1           5           5                              
#> 2           6           6             32             53
#> 3           7           7             75             47
```
