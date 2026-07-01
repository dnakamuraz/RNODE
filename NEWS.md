# RNODE 0.3.0

### New functions

*   `filterSharedTaxa`: given two input matrices, create a new matrix containing shared terminals.
*   `concatenate`: given two input matrices, create a new matrix concatenating partitions.
*   `tnt2nexus`: convert TNT to NEXUS matrices.
*   `nexus2tnt`: convert NEXUS to TNT matrices.

### Improvements

*   `multiSPR`: maxSPR implemented.
*   `filterInvariants`: automatic detection of input format.
*   `filterMissing`: automatic detection of input format.

# RNODE 0.2.0

### New functions

*   `compareBranchLength`: Given two input trees, compare branch lengths of internal edges (shared clades) and terminal edges (shared leaves).
*   `mapBranchLength`: Given one tree without branch lengths and another tree(s) with branch lengths, map the branch lengths from the latter to the former.
*   `mapSupport`: Given one tree with support values and another tree without support values, map the support values from the former to the latter.
*   `filterInvariants`: Given a matrix, delete characters containing only invariants.
*   `splitNoStates`: Given a morphological matrix, split it based on the number of character-states for MK(v) models.

### Improvements

*   `splitOrdFromUnord`: Deleted the parameter "invariant".

### Tutorials

Available for all functions.
