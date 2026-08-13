library(ape)

constraintTree <- function(input_tree, input_constraint, keep_root = TRUE, output_tree = NULL) {
  if (inherits(input_tree, "phylo")) {
    tree <- input_tree
  } else if (is.character(input_tree) && file.exists(input_tree)) {
    tree <- tryCatch({
      ape::read.nexus(input_tree)
    }, error = function(e) {
      ape::read.tree(input_tree)
    })
  } else {
    stop("input_tree must be an object of class 'phylo' or a valid file path.")
  }
  
  taxa <- tree$tip.label
  all_constrained <- unlist(input_constraint)
  if (!all(all_constrained %in% taxa)) {
    stop("Some taxa in input_constraint are not present in the input_tree.")
  }
  
  if (keep_root) {
    root_node <- length(taxa) + 1
    if (root_node <= max(tree$edge)) {
      children <- tree$edge[tree$edge[,1] == root_node, 2]
      
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
    ape::write.nexus(out_phylo, file = output_tree)
  }
  
  return(out_phylo)
}

t <- read.tree(text="((Human,Seal),(Cow,Whale),Mouse);")
res1 <- constraintTree(t, list(c("Human", "Seal")), keep_root=FALSE)
print(write.tree(res1))

res2 <- constraintTree(t, list(c("Human", "Seal", "Cow")), keep_root=FALSE)
print(write.tree(res2))
