# RNODE: Comparisons of topologies, support, and branch lengths between phylogenetic trees

[![language](https://img.shields.io/badge/language-R-blue?style=flat&logo=r&logoColor=white)](https://www.r-project.org)
[![author](https://img.shields.io/badge/author-DYM_Nakamura-blue?logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=c0W8Cm8AAAAJ&hl=en)
[![license](https://img.shields.io/badge/license-GPL_v3-blue?logo=gnu&logoColor=white)](https://www.gnu.org/licenses/gpl-3.0.html)

**RNODE** is an R package to facilitate pre- and postprocessing of phylogenetic analyses, including (1) comparisons of topologies, branch lengths, support values, (2) comparison of DNA sequences, (3) manipulation of cladistic matrices, and (4) manipulation of trees.

Copyright (C) Daniel Y. M. Nakamura 2026

## Cite

If you use **RNODE**, please cite this repository.

## Installation

**RNODE** was tested in R. v. 4.5.2 and can be installed with the following command:

```
devtools::install_github("dnakamuraz/RNODE")
```

## Usage

The following functions are available in **RNODE**:

| Function                  | Class            | Description |
|:--------------------------|:-----------------|:------------|
| *compareBranchLength*     | Tree comparison  | Given two input trees, compare branch lengths of internal edges (shared clades) and terminal edges (shared leaves). The output is a dataframe with node labels and branch lengths.  |
| *multiCID*                | Tree comparison  | Given two sets of binary trees (e.g. MPTs), compute CI distances (normalized or not) between two randomly selected trees or between multiple pairs of trees (all trees or a subsample of them) and then summarize distances as mean, minimum, or maximum values). |
| *multiRF*                 | Tree comparison  | Given two sets of binary trees (e.g. MPTs), compute RF distances (normalized or not) between two randomly selected trees or between multiple pairs of trees (all trees or a subsample of them) and then summarize distances as mean, minimum, or maximum values). |
| *multiSPR*                | Tree comparison  | Given two sets of binary trees (e.g. MPTs), compute SPR distances (normalized or not) between two randomly selected trees or between multiple pairs of trees (all trees or a subsample of them) and then summarize distances as mean, minimum, or maximum values). |
| *normalizedSPR*           | Tree comparison  | Given two binary trees, compute the normalized SPR distance using the upper bound from Ding et al. (2011). |
| *retrodictNodes*          | Tree comparison  | Given two input trees, create a dataframe containing support values of one tree and clade occurrence  of another tree. |
| *sharedNodes*             | Tree comparison  | Given two input trees, compare shared clades. The output is (1) basic statistics about the number of shared clades, support values and their correlation; (2) a dataframe with node labels, descendants, and support values of shared clades, which facilitates descriptive and statistical comparisons of clade composition and support between corresponding nodes.  |
| *summaryTopologicalDist*  | Tree comparison  | Given two sets of trees, compute the number of shared clades, number of unique clades in each tree, Robinson-Foulds, and Cluster Information distance.  |
| *uniqueNodes*             | Tree comparison  | Given two input trees, identify unique clades. The output is two lists containing unique clades and support values in each tree.  |
| *addREP*                  | Tree handling    | Given one alignment and one tree with Goodman-Bremer support values, compute the ratio of explanatory power (REP). |
| *mapBranchLength*         | Tree handling    | Given one tree without branch lengths (e.g. strict consensus) and another tree(s) with branch lengths (e.g. MPTs), map the branch lengths from the latter to the former. |
| *mapSupport*              | Tree handling    | Given one tree with support values (e.g. majority consensus of bootstrap trees) and another tree without support values (e.g. strict consensus of optimal trees), map the support values from the former to the latter. |
| *filterEqualLength*       | Matrix handling  | Given multiple gene alignments, identify gap and gapless files, and write a template of a script considering gap files as unaligned and gapless files as prealigned for POY/PhyG. |
| *concatenate*             | Matrix handling  | Given two molecular or morphological matrices, concatenate them by columns. Taxa present in only one matrix are filled with missing data ('?'). |
| *filterInvariants*        | Matrix handling  | Given a matrix, delete characters containing only invariants. |
| *filterMissing*           | Matrix handling  | Given a matrix, delete taxa and/or characters containing only missing data (?). |
| *filterSharedTaxa*        | Matrix handling  | Given two molecular or morphological matrices, filter taxa using shared taxa only, shared taxa plus taxa unique to input 1, or shared taxa plus taxa unique to input 2. |
| *splitNoStates*           | Matrix handling  | Given a morphological matrix, split it based on the number of character-states for MK(v) models. |
| *splitOrdFromUnord*       | Matrix handling  | Given a morphological matrix and a list of ordered and unordered characters, split the matrix into two matrices. |
| *nexus2tnt*               | Matrix handling  | Converts a NEXUS matrix into TNT format while preserving header data types and ordered character structures. |
| *tnt2nexus*               | Matrix handling  | Converts a TNT matrix into NEXUS format while preserving header data types and ordered character structures. |

The following examples are designed for users with little experience. If you have questions, send a message using GitHub issues. 

#### Tutorials and Examples

For detailed tutorials and examples on how to use these functions, please refer to the package vignettes available on the website's Articles section:

- [Tree Comparisons](https://dnakamuraz.github.io/RNODE/articles/tree-comparisons.html): Identifying shared and unique clades, support comparisons, logistic regressions, branch lengths, and topological distances.
- [Matrix Handling](https://dnakamuraz.github.io/RNODE/articles/matrix-handling.html): Morphological filtering, shared taxa filtering, concatenating matrices, and formatting conversions (TNT / NEXUS).
- [Tree Manipulation](https://dnakamuraz.github.io/RNODE/articles/tree-manipulation.html): Mapping support values and branch lengths across different trees.
