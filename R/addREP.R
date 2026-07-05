#' @title addREP
#' @name addREP
#' @description Given a tree with Goodman-Bremer (GB) support values and the length of the worst tree, compute the ratio of explanatory power (REP).
#' @author Daniel YM Nakamura
#' @param tree A \code{phylo} object that can be loaded using \code{ape::read.tree} for NEWICK files or \code{TreeTools::ReadTntTree} for TNT files.
#' @param best_tree_length The length of the best tree.
#' @param worst_tree_length The length of the worst tree. It can be computed in POY using negative weights.
#' @param output_file A character string specifying the name of the output file in NEXUS format. Default is \code{NULL} (no file is written).
#' @details
#' Grant & Kluge (2007, 2010) proposed REP as a new measure of optimality-based support, in which support values follow the same rank order from GB. REP has the advantage of scaling GB values by the maximum GB, making REP comparable across datasets. REP is calculated as:
#' \deqn{REP = \frac{S' - S}{X - S} = \frac{\mathrm{GB}}{\mathrm{GB}_{\max}}}
#' where S is the optimal length, S' is the length of the tree without a given clade, and X is the length of the worst tree. The maximum value of GB can be calculated using the least parsimonious binary tree (obtained by searching with all characters weighted -1), as this allow groups in the worst tree to contradict with the groups in the best tree.
#' @return A \code{phylo} object with the original Goodman-Bremer support values replaced by REP values.
#' @references Grant, T., & Kluge, A. G. (2007). Ratio of explanatory power (REP): a new measure of group support. Molecular Phylogenetics and Evolution, 44(1), 483-487.
#' @references Grant, T., & Kluge, A. G. (2010). REP provides meaningful measurement of support across datasets. Molecular Phylogenetics and Evolution, 55(1), 340-342.
#' @export
addREP <- function(tree,
                   best_tree_length,
                   worst_tree_length,
                   output_file = NULL) {
    # Check if the object is of class 'phylo'
    if (!inherits(tree, "phylo")) {
        stop("Input 'tree' must be a 'phylo' object.")
    }

    # Check if the tree contains node labels (where GB support values are usually stored)
    if (is.null(tree$node.label)) {
        stop("The input tree does not contain node labels. Ensure GB support values are stored in tree$node.label.")
    }

    # Calculate the denominator (GB_max = X - S)
    GB_max <- worst_tree_length - best_tree_length

    if (GB_max <= 0) {
        stop("worst_tree_length must be strictly greater than best_tree_length.")
    }

    # Safely extract GB support values (convert to numeric, handle NAs/blanks gracefully)
    gb_values <- suppressWarnings(as.numeric(tree$node.label))

    # Calculate the REP values
    rep_values <- gb_values / GB_max

    # Update valid numeric labels with calculated REP values (rounded for readability, e.g., 4 decimal places)
    valid_indices <- !is.na(gb_values)
    tree$node.label[valid_indices] <- as.character(round(rep_values[valid_indices], 4))

    # Write the tree to a NEXUS file if an output_file string is provided
    if (!is.null(output_file)) {
        if (!requireNamespace("ape", quietly = TRUE)) {
            stop("Package 'ape' is required to write NEXUS files. Please install it.")
        }
        ape::write.nexus(tree, file = output_file)
        message(paste("Tree with REP values successfully written to:", output_file))
    }

    # Return the updated tree object
    return(tree)
}
