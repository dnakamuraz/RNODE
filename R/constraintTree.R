#' constraintTree
#'
#' @description Generates a constraint tree based on a list of clades that should be kept monophyletic, collapsing all other clades. This function is useful to generate constraint trees for POY analyses.
#'
#' @param input_tree A phylogenetic tree of class \code{phylo}, or a path to a local NEXUS/Newick file containing the tree.
#' @param input_constraint A list of clades (character vectors of tip labels) that should be kept monophyletic, or a character string of constraints (e.g., \code{"(Taxon1, Taxon2)"}), or a path to a text file containing such constraints. Glob-style wildcards are supported in taxon names: \code{*} matches any characters and \code{?} matches a single character. For example, \code{"(Hylodes*)"} matches all tips whose name starts with \code{"Hylodes"}.
#' @param keep_root A boolean parameter informing whether the root node should be maintained or not even though it does not satisfy \code{input_constraint}.
#' @param output_tree An optional character string specifying the file path for the output tree to be written.
#' @param format A character string specifying the format of the output tree file. Can be either \code{"newick"} (default) or \code{"nexus"}.
#'
#' @return An object of class \code{phylo} representing the constraint tree.
#'
#' @examples
#' library(ape)
#' tree <- read.tree(text = "((Human,Seal),(Cow,Whale),Mouse);")
#'
#' # Using a list of character vectors:
#' constraintTree(tree, input_constraint = list(c("Human", "Seal")))
#'
#' # Using a Newick-style character string for constraints:
#' constraintTree(tree, input_constraint = "(Human,Seal,Cow)")
#'
#' # Using glob wildcards to match multiple taxa:
#' tree2 <- read.tree(text = "((Hylodes1,Hylodes2),(Megaelosia1,Megaelosia2),Outgroup);")
#' constraintTree(tree2, input_constraint = "(Hylodes*,Megaelosia*)")
#'
#' @export
constraintTree <- function(input_tree, input_constraint, keep_root = TRUE, output_tree = NULL, format = c("newick", "nexus")) {
  format <- match.arg(format)

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }

  if (inherits(input_tree, "phylo")) {
    tree <- input_tree
  } else if (is.character(input_tree) && file.exists(input_tree)) {
    tree <- tryCatch(
      {
        ape::read.nexus(input_tree)
      },
      error = function(e) {
        ape::read.tree(input_tree)
      }
    )
  } else {
    stop("input_tree must be an object of class 'phylo' or a valid file path.")
  }

  if (is.character(input_constraint)) {
    if (length(input_constraint) == 1 && file.exists(input_constraint)) {
      input_constraint <- readLines(input_constraint, warn = FALSE)
      input_constraint <- input_constraint[trimws(input_constraint) != ""]
    }

    input_constraint_str <- paste(input_constraint, collapse = ",")

    # Expand glob/regex patterns against tree tip labels
    tokens <- trimws(unlist(strsplit(input_constraint_str, "[(),;]")))
    tokens <- unique(tokens[tokens != ""])
    for (token in tokens) {
      if (grepl("[*?]", token)) {
        # Convert glob wildcards to regex: escape regex metacharacters first, then convert * and ?
        pattern <- gsub("([.+^${}|\\[\\]\\\\])", "\\\\\\1", token)
        pattern <- gsub("\\*", ".*", pattern)
        pattern <- gsub("\\?", ".", pattern)
        pattern <- paste0("^", pattern, "$")
        matches <- grep(pattern, tree$tip.label, value = TRUE)
        if (length(matches) > 0) {
          input_constraint_str <- gsub(token, paste(matches, collapse = ","), input_constraint_str, fixed = TRUE)
        } else {
          warning(paste("No taxa in input_tree matched the pattern:", token))
        }
      }
    }

    text_to_read <- gsub(";;", ";", paste0("(", gsub(";", "", input_constraint_str), ");"))

    parsed_tree <- tryCatch(
      {
        ape::read.tree(text = text_to_read)
      },
      error = function(e) NULL
    )

    if (!is.null(parsed_tree)) {
      clades <- list()
      ntips <- length(parsed_tree$tip.label)
      if (parsed_tree$Nnode > 1) {
        for (i in (ntips + 2):(ntips + parsed_tree$Nnode)) {
          clades[[length(clades) + 1]] <- ape::extract.clade(parsed_tree, i)$tip.label
        }
      } else if (parsed_tree$Nnode == 1) {
        clades[[1]] <- parsed_tree$tip.label
      }
      input_constraint <- clades
    } else {
      input_constraint <- list(trimws(strsplit(gsub("[()]", "", input_constraint_str), ",")[[1]]))
    }
  } else if (!is.list(input_constraint)) {
    input_constraint <- list(input_constraint)
  }

  taxa <- tree$tip.label
  all_constrained <- unlist(input_constraint)
  if (!all(all_constrained %in% taxa)) {
    missing_taxa <- unique(all_constrained[!all_constrained %in% taxa])
    stop(paste("Some taxa in input_constraint are not present in the input_tree:", paste(missing_taxa, collapse = ", ")))
  }

  if (keep_root) {
    root_node <- length(taxa) + 1
    if (root_node <= max(tree$edge)) {
      children <- tree$edge[tree$edge[, 1] == root_node, 2]

      get_tips <- function(tr, node) {
        if (node <= length(tr$tip.label)) {
          return(tr$tip.label[node])
        } else {
          return(ape::extract.clade(tr, node)$tip.label)
        }
      }

      root_clades <- lapply(children, function(child) get_tips(tree, child))
      input_constraint <- c(input_constraint, root_clades)
    }
  }

  # Unique and filter out single-taxon clades
  input_constraint <- unique(lapply(input_constraint, sort))
  input_constraint <- input_constraint[sapply(input_constraint, length) > 1]

  nodes <- as.list(taxa)
  names(nodes) <- taxa
  node_tips <- as.list(taxa)
  names(node_tips) <- taxa

  if (length(input_constraint) > 0) {
    sizes <- sapply(input_constraint, length)
    input_constraint <- input_constraint[order(sizes)]

    for (i in seq_along(input_constraint)) {
      cl <- input_constraint[[i]]
      to_merge <- character()
      for (n in names(node_tips)) {
        if (all(node_tips[[n]] %in% cl)) {
          to_merge <- c(to_merge, n)
        }
      }

      if (length(to_merge) > 1) {
        new_node_name <- paste0("node_", i)
        new_node_str <- paste0("(", paste(unlist(nodes[to_merge]), collapse = ","), ")")
        new_node_tips <- unlist(node_tips[to_merge])

        nodes[to_merge] <- NULL
        node_tips[to_merge] <- NULL

        nodes[[new_node_name]] <- new_node_str
        node_tips[[new_node_name]] <- new_node_tips
      }
    }
  }

  tree_str <- paste0("(", paste(unlist(nodes), collapse = ","), ");")
  out_phylo <- ape::read.tree(text = tree_str)

  if (!is.null(output_tree)) {
    if (format == "nexus") {
      ape::write.nexus(out_phylo, file = output_tree)
    } else {
      ape::write.tree(out_phylo, file = output_tree)
    }
  }

  return(out_phylo)
}
