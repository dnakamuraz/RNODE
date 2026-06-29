# multiSPR

`multiSPR` computes SPR distances between two sets of binary trees (e.g.
MPTs), \\(T_1 = \\{\text{Tree}\_1, \text{Tree}\_2, \dots,
\text{Tree}\_n\\}\\) and \\(T_2 = \\{\text{Tree}\_a, \text{Tree}\_b,
\dots, \text{Tree}\_z\\}\\). The methods available are (1) randomly
selecting one of the binary trees from each set (quick and naive) or (2)
estimating the mean SPR value (or minimum and maximum) from \\(n\\)
pairwise combinations between the two sets. This function is useful when
the two strict consensus trees exhibit polytomies. Both trees must
contain the same set of leaves.

## Usage

``` r
multiSPR(
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

  Optional. Specify if SPR distances will be calculated by (1) `random`
  (default: selects one binary tree randomly from the multiPhylo
  object), (2) `meanSPR` (calculates mean of all pairwise SPR distances
  between two `multiPhylo` objects), (3) `minSPR` (calculates the
  minimum value), (4) `maxSPR` (calculates the maximum value), and (5)
  `all` (calculates the mean, minimum, and maximum values).

- normalization:

  Optional. Specify if SPR distances should be normalized using upper
  bound values (Ding et al. 2011). See details in
  [normalizedSPR](https://dnakamuraz.github.io/RNODE/reference/normalizedSPR.md).
  By default, SPR distances are not normalized.

- subsample:

  Optional. Specify if SPR distances will be calculated to only a
  fraction of the total number of trees available in each set (default =
  1; i.e. all trees are evaluated). Zero values are not accepted.

## References

Ding, Y., Grünewald, S., Humphries, P.J., 2011. On agreement forests. J.
Comb. Theory Ser. 118(7), 2059–2065.

## Author

Daniel YM Nakamura
