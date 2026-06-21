#' @title addREP
#' @name addREP
#' @description Given one alignment and one tree with Goodman-Bremer (GB) support values, compute the ratio of explanatory power (REP) using equally weighted parsimony.
#' @author Daniel YM Nakamura
#' @param data A cladistic matrix.
#' @param tree A \code{phylo} object that can be loaded using \code{ape::read.tree} for NEWICK files or \code{TreeTools::ReadTntTree} for TNT files.
#' @details
#' Grant & Kluge (2007, 2010) proposed REP as a new measure of optimality-based support, in which support values follow the same rank order from GB. REP has the advantage of scaling GB values by the maximum GB, making REP comparable across datasets. REP is calculated as:
#' \deqn{REP = \frac{S' - S}{X - S} = \frac{\mathrm{GB}}{\mathrm{GB}_{\max}}}
#' where S is the optimal length, S' is the length of the tree without a given clade, and X is the length of the worst tree. The maximum value of GB can be calculated using the least parsimonious binary tree (obtained by searching with all characters weighted -1), as this allow groups in the worst tree to contradict with the groups in the best tree.
#' @references Grant, T., & Kluge, A. G. (2007). Ratio of explanatory power (REP): a new measure of group support. Molecular Phylogenetics and Evolution, 44(1), 483-487.
#' @references Grant, T., & Kluge, A. G. (2010). REP provides meaningful measurement of support across datasets. Molecular Phylogenetics and Evolution, 55(1), 340-342.
#' @export
addREP = function(data,
               tree){

}


#data("Laurasiatherian")
#alignment = Laurasiatherian
#dna_states <- c("a", "c", "g", "t")
#cost_matrix <- matrix(
#    +     c(0,-1,-1,-1,
#            +       -1,0,-1,-1,
#            +       -1,-1,0,-1,
#            +       -1,-1,-1,0),
#    +     nrow = 4, byrow = TRUE,
#    +     dimnames = list(dna_states, dna_states)
#    + )
#pratchet(alignment, cost=cost_matrix)
