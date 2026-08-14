# constraintTree

Generates a constraint tree based on a list of clades that should be
kept monophyletic, collapsing all other clades.

## Usage

``` r
constraintTree(
  input_tree,
  input_constraint,
  keep_root = TRUE,
  output_tree = NULL,
  format = c("nexus", "newick")
)
```

## Arguments

- input_tree:

  A phylogenetic tree of class `phylo`, or a path to a local
  NEXUS/Newick file containing the tree.

- input_constraint:

  A list of clades (character vectors of tip labels) that should be kept
  monophyletic, or a character string of constraints (e.g., "(Taxon1,
  Taxon2)"), or a path to a text file containing such constraints.

- keep_root:

  A boolean parameter informing whether the root node should be
  maintained or not even though it does not satisfy `input_constraint`.

- output_tree:

  An optional character string specifying the file path for the output
  tree to be written.

- format:

  A character string specifying the format of the output tree file. Can
  be either `"nexus"` (default) or `"newick"`.

## Value

An object of class `phylo` representing the constraint tree.

## Examples

``` r
library(ape)
tree <- read.tree(text = "((Human,Seal),(Cow,Whale),Mouse);")

# Using a list of character vectors:
constraintTree(tree, input_constraint = list(c("Human", "Seal")))
#> 
#> Phylogenetic tree with 5 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Mouse, Human, Seal, Cow, Whale
#> 
#> Unrooted; no branch length.

# Using a Newick-style character string for constraints:
constraintTree(tree, input_constraint = "(Human,Seal,Cow)")
#> 
#> Phylogenetic tree with 5 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Mouse, Human, Seal, Cow, Whale
#> 
#> Unrooted; no branch length.
```
