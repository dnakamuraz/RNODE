# multiRF

`multiRF` computes Robinson-Foulds distances between two sets of binary
trees (e.g. MPTs), \\(T_1 = \\{\text{Tree}\_1, \text{Tree}\_2, \dots,
\text{Tree}\_n\\}\\) and \\(T_2 = \\{\text{Tree}\_a, \text{Tree}\_b,
\dots, \text{Tree}\_z\\}\\). The methods available are (1) randomly
selecting one of the binary trees from each set (quick and naive) and
(2) estimating the mean RF (or minimum or maximum) from \\(n\\) pairwise
combinations between the two sets. Both trees must contain the same set
of leaves.

## Usage

``` r
multiRF(
  trees1,
  trees2,
  method = "random",
  normalization = FALSE,
  subsample = 1
)
```

## Arguments

- trees1:

  A `phylo` or `multiPhylo` object with multiple trees that can be
  loaded using
  [`ape::read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html) for
  NEWICK files or
  [`TreeTools::ReadTntTree`](https://ms609.github.io/TreeTools/reference/ReadTntTree.html)
  for TNT files. If the pool of MPTs presents binary and non-binary
  trees, only binary trees are processed.

- trees2:

  Another `phylo` or `multiPhylo` object.

- method:

  Optional. Specify if RF distances will be calculated by (1) `random`
  (default: selects one binary tree randomly from the multiPhylo
  object), (2) `meanRF` (calculates mean of all pairwise RF distances
  between two `multiPhylo` objects), (3) `minRF` (calculates the
  minimum), (4) `maxRF` (calculates the maximum), or (5) `all`
  (calculates mean, minimum, and maximum values).

- normalization:

  Optional. Specify if RF distances should be normalized using the
  maximum possible RF distance \\(i(T_1) + i(T_2)\\) (sum of clades in
  both trees). By default, RF distances are not normalized.

- subsample:

  Optional. Specify if RF distances will be calculated to only a
  fraction of the total number of trees available in each set (default =
  1; i.e. all trees are evaluated). Zero values are not accepted.

## Author

Daniel YM Nakamura
