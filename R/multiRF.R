#' @title multiRF
#' @name multiRF
#' @description \code{multiRF} computes Robinson-Foulds distances between two sets of binary trees (e.g. MPTs), \eqn{T_1 = \{\text{Tree}_1, \text{Tree}_2, \dots, \text{Tree}_n\}} and \eqn{T_2 = \{\text{Tree}_a, \text{Tree}_b, \dots, \text{Tree}_z\}}. The methods available are (1) randomly selecting one of the binary trees from each set (quick and naive) and (2) estimating the mean RF (or minimum or maximum) from \eqn{n} pairwise combinations between the two sets. Both trees must contain the same set of leaves.
#' @author Daniel YM Nakamura
#'
#' @param trees1 A \code{phylo} or \code{multiPhylo} object with multiple trees that can be loaded using \code{ape::read.tree} for NEWICK files or \code{TreeTools::ReadTntTree} for TNT files. If the pool of MPTs presents binary and non-binary trees, only binary trees are processed.
#' @param trees2 Another \code{phylo} or \code{multiPhylo} object.
#' @param method Optional. Specify if RF distances will be calculated by (1) \code{random} (default: selects one binary tree randomly from the multiPhylo object), (2) \code{meanRF} (calculates mean of all pairwise RF distances between two \code{multiPhylo} objects), (3) \code{minRF} (calculates the minimum), (4) \code{maxRF} (calculates the maximum), or (5) \code{all} (calculates mean, minimum, and maximum values).
#' @param normalization Optional. Specify if RF distances should be normalized using the maximum possible RF distance \eqn{i(T_1) + i(T_2)} (sum of clades in both trees). By default, RF distances are not normalized.
#' @param subsample Optional. Specify if RF distances will be calculated to only a fraction of the total number of trees available in each set (default = 1; i.e. all trees are evaluated). Zero values are not accepted.
#'
#' @export
multiRF = function(trees1, trees2,
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

    rf_distances <- matrix(NA, nrow = length(trees1), ncol = length(trees2))

    for (i in seq_along(trees1)) {
      for (j in seq_along(trees2)) {

        rf_distances[i, j] <- phangorn::RF.dist(
          trees1[[i]], trees2[[j]],
          normalize = normalization
        )

      }
    }

    vals <- as.vector(rf_distances)

    if (method == "meanRF") {
      return(mean(vals, na.rm = TRUE))

    } else if (method == "minRF") {
      return(min(vals, na.rm = TRUE))

    } else if (method == "maxRF") {
      return(max(vals, na.rm = TRUE))

    } else if (method == "random") {
      return(sample(vals, 1))

    } else if (method == "all") {
      return(data.frame(
        minRF = min(vals, na.rm = TRUE),
        maxRF = max(vals, na.rm = TRUE),
        meanRF = mean(vals, na.rm = TRUE)
      ))
    }
  }

  ##########
  # CASE 2 #
  ##########
  else if (inherits(trees1, "multiPhylo") && inherits(trees2, "phylo")) {

    distances <- sapply(trees1, function(tree) {
      phangorn::RF.dist(tree, trees2, normalize = normalization)
    })

    if (method == "meanRF") return(mean(distances))
    else if (method == "minRF") return(min(distances))
    else if (method == "maxRF") return(max(distances))
    else if (method == "random") return(sample(distances, 1))
    else if (method == "all") {
      return(data.frame(
        minRF = min(distances),
        maxRF = max(distances),
        meanRF = mean(distances)
      ))
    }
  }

  ##########
  # CASE 3 #
  ##########
  else if (inherits(trees1, "phylo") && inherits(trees2, "multiPhylo")) {

    distances <- sapply(trees2, function(tree) {
      phangorn::RF.dist(trees1, tree, normalize = normalization)
    })

    if (method == "meanRF") return(mean(distances))
    else if (method == "minRF") return(min(distances))
    else if (method == "maxRF") return(max(distances))
    else if (method == "random") return(sample(distances, 1))
    else if (method == "all") {
      return(data.frame(
        minRF = min(distances),
        maxRF = max(distances),
        meanRF = mean(distances)
      ))
    }
  }

  ##########
  # CASE 4 #
  ##########
  else if (inherits(trees1, "phylo") && inherits(trees2, "phylo")) {

    return(phangorn::RF.dist(trees1, trees2, normalize = normalization))
  }

  stop("Invalid tree input types")
}
