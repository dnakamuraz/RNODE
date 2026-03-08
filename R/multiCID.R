#' @title multiCID
#' @name multiCID
#' @description \code{multiCID} computes CID distances between two sets of binary trees (e.g. MPTs), \eqn{T_1 = \{\text{Tree}_1, \text{Tree}_2, \dots, \text{Tree}_n\}} and \eqn{T_2 = \{\text{Tree}_a, \text{Tree}_b, \dots, \text{Tree}_z\}}. The methods available are (1) randomly selecting one of the binary trees from each set (quick and naive) and (2) estimating the mean CID (or minimum or maximum) from \eqn{n} pairwise combinations between the two sets. Both trees must contain the same set of leaves.
#' @author Daniel YM Nakamura
#'
#' @param trees1 A \code{phylo} or \code{multiPhylo} object with multiple trees that can be loaded using \code{ape::read.tree} for NEWICK files or \code{TreeTools::ReadTntTree} for TNT files. If the pool of MPTs presents binary and non-binary trees, only binary trees are processed.
#' @param trees2 Another \code{phylo} or \code{multiPhylo} object.
#' @param method Optional. Specify if CID distances will be calculated by (1) \code{random} (default: selects one binary tree randomly from the multiPhylo object), (2) \code{meanCID} (calculates mean of all pairwise CID distances between two \code{multiPhylo} objects), (3) \code{minCID} (calculates the minimum), (4) \code{maxCID} (calculates the maximum), or (5) \code{all} (calculates mean, minimum, and maximum values).
#' @param normalization Optional. Specify if CID distances should be normalized. By default, CID distances are not normalized.
#' @param subsample Optional. Specify if CID distances will be calculated to only a fraction of the total number of trees available in each set (default = 1; i.e. all trees are evaluated). Zero values are not accepted.
#'
#'
#' @export
multiCID = function(trees1, trees2,
                    method = "random",
                    normalization = FALSE,
                    subsample = 1) {

  if (subsample <= 0) stop("subsample must be > 0")
  if (subsample > 1) stop("subsample must be <= 1")

  # Keep only binary trees
  if (inherits(trees1, "multiPhylo")) {
    trees1 <- trees1[sapply(trees1, is.binary)]
  }
  if (inherits(trees2, "multiPhylo")) {
    trees2 <- trees2[sapply(trees2, is.binary)]
  }

  # Subsample trees
  if (subsample < 1) {
    if (inherits(trees1, "multiPhylo")) {
      n1 <- ceiling(length(trees1) * subsample)
      trees1 <- trees1[sample(seq_along(trees1), n1)]
    }
    if (inherits(trees2, "multiPhylo")) {
      n2 <- ceiling(length(trees2) * subsample)
      trees2 <- trees2[sample(seq_along(trees2), n2)]
    }
  }

  ##########
  # CASE 1 #
  ##########
  if (inherits(trees1, "multiPhylo") && inherits(trees2, "multiPhylo")) {

    CID_distances <- matrix(NA, nrow = length(trees1), ncol = length(trees2))

    for (i in seq_along(trees1)) {
      for (j in seq_along(trees2)) {

        CID_distances[i, j] <- ClusteringInfoDistance(
          trees1[[i]], trees2[[j]],
          normalize = normalization
        )

      }
    }

    vals <- as.vector(CID_distances)

    if (method == "meanCID") {
      return(mean(vals, na.rm = TRUE))

    } else if (method == "minCID") {
      return(min(vals, na.rm = TRUE))

    } else if (method == "maxCID") {
      return(max(vals, na.rm = TRUE))

    } else if (method == "random") {
      return(sample(vals, 1))

    } else if (method == "all") {
      return(data.frame(
        minCID = min(vals, na.rm = TRUE),
        maxCID = max(vals, na.rm = TRUE),
        meanCID = mean(vals, na.rm = TRUE)
      ))
    }
  }

  ##########
  # CASE 2 #
  ##########
  else if (inherits(trees1, "multiPhylo") && inherits(trees2, "phylo")) {

    distances <- sapply(trees1, function(tree) {
      ClusteringInfoDistance(tree, trees2, normalize = normalization)
    })

    if (method == "meanCID") return(mean(distances))
    else if (method == "minCID") return(min(distances))
    else if (method == "maxCID") return(max(distances))
    else if (method == "random") return(sample(distances, 1))
    else if (method == "all") {
      return(data.frame(
        minCID = min(distances),
        maxCID = max(distances),
        meanCID = mean(distances)
      ))
    }
  }

  ##########
  # CASE 3 #
  ##########
  else if (inherits(trees1, "phylo") && inherits(trees2, "multiPhylo")) {

    distances <- sapply(trees2, function(tree) {
      ClusteringInfoDistance(trees1, tree,normalize = normalization)
    })

    if (method == "meanCID") return(mean(distances))
    else if (method == "minCID") return(min(distances))
    else if (method == "maxCID") return(max(distances))
    else if (method == "random") return(sample(distances, 1))
    else if (method == "all") {
      return(data.frame(
        minCID = min(distances),
        maxCID = max(distances),
        meanCID = mean(distances)
      ))
    }
  }

  ##########
  # CASE 4 #
  ##########
  else if (inherits(trees1, "phylo") && inherits(trees2, "phylo")) {
    return (ClusteringInfoDistance(trees1, trees2,normalize = normalization))
  }

  stop("Invalid tree input types")
}

'
multiCID = function(trees1, trees2,
                    method="random",
                    normalization=F){


  if (class(trees1)=="multiPhylo"){
    trees1 = trees1[sapply(trees1, is.binary)]
  }
  if (class(trees2)=="multiPhylo"){
    trees2 = trees2[sapply(trees2, is.binary)]
  }


  #################
  # CID DISTANCES #
  #################

  # 1. If both input files are multiPhylo
  if (class(trees1)=="multiPhylo" && class(trees2)=="multiPhylo"){
    # 1.1 If the method is "meanCID"
    if (method=="meanCID"){
      # Create empty matrix
      CID_distances <- matrix(NA, nrow = length(trees1), ncol = length(trees2))
      # Compute pairwise CID distances
      for (i in seq_along(trees1)){
        for (j in seq_along(trees2)){
          # 1.1.1 If normalization is required
          if (normalization==T){CID_distances[i, j] = RNODE::normalizedCID(trees1[[i]], trees2[[j]])}
          # 1.1.2 Else, if normalization is not required
          else if (normalization==F){CID_distances[i, j] = TreeDist::CIDDist(trees1[[i]], trees2[[j]])}
        }
      }
      # Calculate mean CID from all values in the matrix
      meanCID <- mean(CID_distances, na.rm = TRUE)
      return(list(meanCID, CID_distances))
    }

    # 1.2 If the method is "minCID"
    if (method=="minCID"){
      # Create empty matrix
      CID_distances <- matrix(NA, nrow = length(trees1), ncol = length(trees2))
      # Compute pairwise CID distances
      for (i in seq_along(trees1)){
        for (j in seq_along(trees2)){
          # 1.1.1 If normalization is required
          if (normalization==T){CID_distances[i, j] = RNODE::normalizedCID(trees1[[i]], trees2[[j]])}
          # 1.1.2 Else, if normalization is not required
          else if (normalization==F){CID_distances[i, j] = TreeDist::CIDDist(trees1[[i]], trees2[[j]])}
        }
      }
      # Calculate mean CID from all values in the matrix
      minCID <- min(CID_distances, na.rm = TRUE)
      return(list(minCID, CID_distances))
    }

    # 1.3 If the method is "random"
    else if (method=="random"){
      # Randomly select one tree 1
      random_tree1 = trees1[[sample(1:length(trees1),1)]]
      # Randomly select one tree 2
      random_tree2 = trees2[[sample(1:length(trees2),1)]]
      # 1.2.1 If normalization is required, compute normalized CID distance
      if (normalization==T) {result = RNODE::normalizedCID(random_tree1, random_tree2)}
      # 1.2.2 If normalization is not required, compute CID distance
      else if (normalization==F) {result = TreeDist::CIDDist(random_tree1,random_tree2)}
      return(result)
    }
  }

  # 2. Else, if one or both of the trees are not multiPhylo
  else {
    # 2.1 If normalization is required
    if (normalization==T){
      result = RNODE::normalizedCID(trees1, trees2)
      if (method=="meanCID"){result = mean(result)}
      else if (method=="minCID"){result = min(result)}
      else if (method=="random"){result = sample(result,1)}
      return(result)
    }
    # 2.2 If normalization is not required
    else if (normalization==F){
      result = TreeDist::CIDDist(trees1, trees2)
      if (method=="meanCID"){result = mean(result)}
      else if (method=="minCID"){result = min(result)}
      else if (method=="random"){result = sample(result,1)}
      return(result)
    }
  }
}


#TESTE
trees1 = mp_mol_mpts[[2]]$mp_mol_mpts
trees2 = mp_te_mpts[[2]]$mp_te_mpts
trees1
trees2

shared_terminals <- intersect(trees1$tip.label, trees2$tip.label)
setdiff(trees1$tip.label, trees2$tip.label) # unique terminals in tree 1
setdiff(trees2$tip.label, trees1$tip.label) # unique terminals in tree 2
trees1 <- drop.tip(trees1, trees1$tip.label[!(trees1$tip.label %in% shared_terminals)])
trees2 <- drop.tip(trees2, trees2$tip.label[!(trees2$tip.label %in% shared_terminals)])

CID.dist(trees1, trees2)
normalizedCID(trees1, trees2)
multiCID(mp_mol_mpts[[24]]$mp_mol_mpts, mp_te_mpts[[24]]$mp_te_mpts, method="minCID", normalization=T)

a = vector("list", 1) # empty list
for (i in 7) {a[[i]] = multiCID(ml_mol_best[[i]],
                       ml_te_noasc_best[[i]],
                       method="meanCID",
                       normalization=T)}
a

'
