# addREP

Given one alignment and one tree with Goodman-Bremer (GB) support
values, compute the ratio of explanatory power (REP) using equally
weighted parsimony.

## Usage

``` r
addREP(data, tree)
```

## Arguments

- data:

  A cladistic matrix.

- tree:

  A `phylo` object that can be loaded using
  [`ape::read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html) for
  NEWICK files or
  [`TreeTools::ReadTntTree`](https://ms609.github.io/TreeTools/reference/ReadTntTree.html)
  for TNT files.

## Details

Grant & Kluge (2007, 2010) proposed REP as a new measure of
optimality-based support, in which support values follow the same rank
order from GB. REP has the advantage of scaling GB values by the maximum
GB, making REP comparable across datasets. REP is calculated as: \$\$REP
= \frac{S' - S}{X - S} = \frac{\mathrm{GB}}{\mathrm{GB}\_{\max}}\$\$
where S is the optimal length, S' is the length of the tree without a
given clade, and X is the length of the worst tree. The maximum value of
GB can be calculated using the least parsimonious binary tree (obtained
by searching with all characters weighted -1), as this allow groups in
the worst tree to contradict with the groups in the best tree.

## References

Grant, T., & Kluge, A. G. (2007). Ratio of explanatory power (REP): a
new measure of group support. Molecular Phylogenetics and Evolution,
44(1), 483-487.

Grant, T., & Kluge, A. G. (2010). REP provides meaningful measurement of
support across datasets. Molecular Phylogenetics and Evolution, 55(1),
340-342.

## Author

Daniel YM Nakamura
