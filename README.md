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
| *mapBranchLength*         | Tree handling    | Given one tree without branch lengths (e.g. strict consensus) and another tree(s) with branch lengths (e.g. MPTs), map the branch lengths from the latter to the former. |
| *mapSupport*              | Tree handling    | Given one tree with support values (e.g. majority consensus of bootstrap trees) and another tree without support values (e.g. strict consensus of optimal trees), map the support values from the former to the latter. |
| *rep*                     | Tree handling    | Given one alignment and one tree with Goodman-Bremer support values, compute the ratio of explanatory power (REP). |
| *findEqualLength*         | Matrix handling  | Given multiple gene alignments, identify gap and gapless files, and write a template of a script considering gap files as unaligned and gapless files as prealigned for POY/PhyG. |
| *filterInvariants*        | Matrix handling  | Given a matrix, delete characters containing only invariants. |
| *filterMissing*           | Matrix handling  | Given a matrix, delete taxa and/or characters containing only missing data (?). |
| *splitNoStates*           | Matrix handling  | Given a morphological matrix, split it based on the number of character-states for MK(v) models. |
| *splitOrdFromUnord*       | Matrix handling  | Given a morphological matrix and a list of ordered and unordered characters, split the matrix into two matrices. |

The following examples are designed for users with little experience. If you have questions, send a message using GitHub issues. 

### Example 1: Tree comparisons

#### Example 1.1 Identify shared and unique clades

Using simple simulations, we can demonstrate how to compare support values between trees. We first simulate two trees containing support values: 

```
set.seed(44)
# Simulate tree a
a = pbtree(n=7) 
# Generate random support values as integers to tree a
node_labels = sample(1:100, a$Nnode, replace = TRUE) 
# Add the support values as node labels to tree a
a$node.label = node_labels 

set.seed(88)
# Simulate tree b
b = pbtree(n=7)
# Generate random support values as integers to tree b
node_labels = sample(1:100, b$Nnode, replace = TRUE) 
# Add the support values as node labels to tree b
b$node.label = node_labels 
```

Next, we run *sharedNodes* to identify matching clades and their descendants and support values. Additionally, we also can plot the trees.

```
# Compare shared clades and support values (and plot)
df = sharedNodes(tree1=a, tree2=b, composition=T, 
                 plotTrees = T,
                 output.tree1="example1.1_simulated1.pdf",
                 output.tree2="example1.1_simulated2.pdf", 
                 tree.width = 3, # adjust tree width
                 tree.height = 4, # adjust tree height
                 tree.fsize = 1, # adjust font size
                 tree.adj=c(1.2,3), # adjust support position
                 tree.cex=.5, # adjust support size
                 node.numbers=T) # show node index
```

<p align="center">
  <a href="tutorial/example1.1_df.png"><img src="tutorial/example1.1_df.png" alt="df" width="100%"></a>
</p>


<p align="center">
  <a href="tutorial/example1.1_simulated1.png"><img src="tutorial/example1.1_simulated1.png" alt="Fig 1" width="45%"></a>
  <a href="tutorial/example1.1_simulated2.png"><img src="tutorial/example1.1_simulated2.png" alt="Fig 2" width="45%"></a>
</p>

Alternatively, we can identify and plot the unique clades:

```
uniqueNodes(a, b, composition=T, dataframe=T,
            plotTrees=T, output.tree = "example1.1_unique.pdf",
            node.numbers=T, 
            tree.fsize=2, # adjust text size
            tree.cex=7.5, # adjust circle size
            sup.adj1=c(-.2,4), # adjust support from tree a
            sup.adj2=c(1.3,4) # adjust support from tree b
            )
```

<p align="center">
  <a href="tutorial/example1.1_unique.png"><img src="tutorial/example1.1_unique.png" alt="df" width="100%"></a>
</p>

#### Example 1.2 Support comparisons

We can use *sharedNodes* to compare two empirical trees in .nwk format estimated in TNT. Polytomies and input trees with different taxon samples are accepted but names of corresponding leaves should be equal in the input trees. For instance, using the data set from Whitcher et al. (2025), we can plot the relationship of bootstrap values between molecular (MOL) and total evidence (TE) trees analyzed in TNT.

```
# Load trees
MOL = read.tree("../testdata/051b_MOL_BS_TNT.nwk")
TE = read.tree("../testdata/051d_TE_BS_TNT.nwk")

# Run sharedNodes
df = sharedNodes(tree1=MOL, tree2=TE, spearman = T)

# Plot the relationship of support between trees
ggplot(df, aes(as.numeric(Support_Tree_1), as.numeric(Support_Tree_2))) +
  geom_point(size = 5, show.legend = F, alpha=.5) +
  theme_minimal() + 
  geom_smooth(method = "lm", se = T, color = "red", linewidth = .5) +
  labs(x="\n Bootstrap in the MOL tree",
       y="Bootstrap in the TE tree \n")
```

<p align="center">
  <a href="tutorial/example1.2_correlation.png"><img src="tutorial/example1.2_correlation.png" alt="df" width="100%"></a>
</p>

As expected, there is a significant correlation between bootstrap values of MOL and TE trees (Spearman: rho = 0.89; P < 0.001). 

#### Example 1.3 Logistic regressions

If the user wants to test if support values of one tree predict the occurrence of clades in another tree, the function *retrodictNodes* creates a dataframe containing support values of tree 1 and the occurrence of the clade in tree 2, which can be used for logistic regressions.

```
# Load trees
MOL = read.tree("../testdata/001_MOL_IQTREE.contree")
TE = read.tree("../testdata/001_TE_ASC_IQTREE.contree")

# Run retrodictNodes
df = retrodictNodes(MOL, TE)
df$occurrence_tree2 = as.factor(df$occurrence_tree2)

# Fit the logistic regression
model <- glm(occurrence_tree2 ~ support_tree1, data = df, family = binomial)
summary(model)

# Convert log-odds to odd ratios
exp(coef(model))
```

Using the data set from Janssens et al. (2018), the logistic regression revealed an intercept of 0.009 (i.e. when bootstrap is 0 in the first tree, the odds of presence of the clade in the second tree is 0.009; P < 0.01). Furthermore, for every one-unit increase in bootstrap in the first tree, the odds of presence of the clade in the second tree increase by 1.075 (7.5%). 

#### Example 1.4 Branch length comparisons

In addition to descendants and support values, branch lengths can be compared. 

```
# Compare branch lengths
RNODE::compareBranchLength(MOL, TE, composition=T)

# Correlation between branch lengths
summary(lm(data=df, formula=EdgeLength_tree1 ~ EdgeLength_tree2))

# Plot 
ggplot(df, aes(as.numeric(EdgeLength_tree1), as.numeric(EdgeLength_tree2))) +
  geom_point(size = 5, show.legend = F, alpha=.5) +
  theme_minimal() + 
  geom_smooth(method = "lm", se = T, color = "red", linewidth = .5) +
  labs(x="\n Branch lengths in the MOL tree",
       y="Branch lengths in the TE tree \n")
```

<p align="center">
  <a href="tutorial/example1.4_lengths.png"><img src="tutorial/example1.4_lengths.png" alt="df" width="100%"></a>
</p>

As expected, there is a significant correlation between bootstrap values of MOL and TE trees (lsinear model: estimate = 0.78; R-squared = 0.86; P < 0.001).

#### Example 1.5 Topological distances

In addition to comparisons between shared clades, support values and branch lengths, a popular method to compare phylogenies is based on topological metrics. Popular metrics like Robinson-Foulds and Cluster Information distances can be summarized using *summaryTopologicalDist*. Moreover, a common topological metric is the number of SPR moves to edit one tree into another tree. However, implementations are lacking in R to normalize SPR distances using the refined upper bound from Ding et al. (2011) (*normalizedSPR*) and computing SPR distances for multiple trees (*multiSPR*). 

```
# Read trees
mol = read.tree("../testdata/003_MOL_IQTREE.contree")
te = read.tree("../testdata/003_TE_ASC_IQTREE.contree")

# RF and CID
summaryTopologicalDist(mol, te)

# Normalized SPR
normalizedSPR(mol, te)
```

The normalized SPR is 0.1931574.

### Example 2 Comparison of DNA sequences

### Example 3 Matrix handling

#### Example 3.1 Morphological matrix

The function *filterMissing* deletes taxa and/or characters containing only missing data. In the following example, the output file will be saved as *test_filterMissing_FILTERED.nexus*:

```
filterMissing(input="../testdata/test_filterMissing.nexus", 
              input_format="nexus",
              output_path="../testdata/test_filterMissing",
              missing="both")
```

The function *filterInvariants* deletes invariant characters, which is useful to accelerate the graph searches. In Maximum Likelihood and Bayesian analyses using the MKv model with ascertainment bias correction (ASC), invariants must be deleted. Here, we follow the definition of invariant from IQ-Tree, characterized by: (1) constant sites containing only a single character state in all sequences, (2) partially constant sites (N and/or -), and (3) ambiguously constant sites (e.g. C, Y and -). In the following example, 122 invariants are detected.

```
filterInvariants(input="../testdata/015_MORPH_data.nexus",
                 input_format = "nexus",
                 output_index="../testdata/015_MORPH_data")
```

The function *splitOrdFromUnord* splits a morphological matrix into partitions of ordered and unordered characters based on a list of ordered characters.

```
# Data input of list of ordered characters
list_ordered=c(1, 6, 7, 8, 10, 12, 13, 14, 17, 19, 23, 26, 31, 35, 41, 44, 45, 48, 51, 54, 55, 68, 71, 72, 92, 94, 96, 102, 105, 108, 109, 128, 129, 130, 131, 132, 135, 142, 144, 152, 153, 193)

splitOrdFromUnord(input="../testdata/048_MORPH_data.nex", 
                  input_format = "nexus",
                  output_index = "../testdata/048_MORPH", 
                  list_ordered=list_ordered)
```

The function *splitNoStates* splits characters from a morphological matrix according to their number of character-states. This procedure has been recommended to run phylogenetic analyses with the MK and MKv models (the 'K' refers to the number of states). Khakurel et al. (2024) demonstrated that MK models with high K values can understimate the branch lengths, whereas MK models with small K values can overstimate them. As such, some recent studies have partitioned morphological characters according to their number of states (e.g. Černý & Simonoff 2023).

```
plitNoStates(input = "../testdata/015_MORPH_data.nexus", 
             input_format = "nexus", 
             output_index = "../testdata/015_MORPH_data", 
             ambiguity_addState = T, 
             inapplicable_addState = T, 
             log=T, 
             write=T)
```

#### Example 3.2 Dynamic homology



### Example 4 Tree manipulation

#### Example 4.1 Mapping support

Given a tree A without support values (e.g. strict consensus of optimal trees) and a tree B with support values (e.g. majority consensus from bootstrap pseudo-replicates), *mapSupport* returns the tree A with support values from shared clades with tree B. For instance, the strict consensus of optimal trees and the majority consensus tree from bootstrap trees share 223 clades, presenting 6 unique clades in the strict consensus and 1 unique clade in the bootstrap tree.  

```
# Read trees
opt = read.tree("../testdata/051a_strictConsensus_MOL_TNT_results.nwk")
BS = read.tree("../testdata/051b_MOL_BS_TNT.nwk")

# Compute topological distances
summaryTopologicalDist(opt, BS)

# Map the BS values from the majority consensus tree to the optimal tree
opt_with_bs = mapSupport(opt, BS)
opt_with_bs[1]
```

#### Example 4.2 Mapping branch lengths

Finally, we can map branch lengths to the strict consensus either using the minimum values from a pool of MPTs or randomly selecting one of the MPTs. Using empirical data from Nakamura et al. 2025, we demonstrate the function *mapBranchLength*:

```
strict = read.tree("../testdata/cymb_IP_GB.1.nwk")
mpts = read.tree("../testdata/cymb_IP_trees.nwk")

# Map using minimum branch length per shared edge
mapped_min <- mapBranchLength(strict, mpts, method = "minimum")

# Plot strict + mapped side by side
png("../tutorial/example4.2.png", width = 2400, height = 4000, res = 150)
par(mfrow = c(1, 2),            # 1 row, 2 columns
    mar = c(4, 4, 2, 2),        # margins
    oma = c(1, 1, 1, 1))        # outer margins
# Panel 1 — strict consensus
plot(ladderize(strict),
     main = "Strict Consensus Tree",
     cex = 0.8)
# Panel 2 — mapped branch lengths
plot(ladderize(mapped_min),
     main = "Mapped Branch Lengths (Minimum)",
     cex = 0.8)
dev.off()
```

<p align="center">
  <a href="tutorial/example4.2.png"><img src="tutorial/example4.2.png" alt="df" width="100%"></a>
</p>
