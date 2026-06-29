# mapBranchLength

Given a tree 1 without branch lengths (e.g. strict consensus of optimal
trees) and a tree 2 with branch lengths (e.g. each MPT before
reconciliation into a strict consensus), return tree 1 with branch
lengths from shared clades with tree 2.

## Usage

``` r
mapBranchLength(tree1, trees2, method = "minimum")
```

## Arguments

- tree1:

  A phylo object without branch lengths

- trees2:

  A phylo object with branch lengths

- method.:

  Optional. Method to sample branch lengths. If "random", randomly
  select one of the trees from trees2. If "minimum" (default), map the
  minimum values from the pool of MPTs (trees2) to each edge present in
  tree1

## Author

Daniel YM Nakamura

## Examples

``` r
# Example 1 (identify unique nodes)
tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
mapBranchLength (tree1, tree2)
#> 
#> Phylogenetic tree with 5 tips and 4 internal nodes.
#> 
#> Tip labels:
#>   t1, t3, t2, t4, t5
#> 
#> Rooted; includes branch length(s).
```
