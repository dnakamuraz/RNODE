#' @title splitNoStates
#' @name splitNoStates
#' @description Splits a morphological matrix according to the number of character-states. This procedure has been recommended to run phylogenetic analyses with the MK and MKv models (the 'K' refers to the number of states). Khakurel et al. (2024) demonstrated that MK models with high K values can underestimate the branch lengths, whereas MK models with small K values can overstimate them. As such, some recent studies have partitioned morphological characters according to their number of states (e.g. Černý & Simonoff 2023).
#' @author Daniel YM Nakamura
#' @param input Input file (morphological matrix, either loaded previously in R or loaded from a local file).
#' @param inpu_format To load from a local file: 'nexus' or 'tnt'.
#' @param ambiguity_addState If T, ambiguities are counted as additional character-states (default: ambiguity_addState = F).
#' @param inapplicable_addState If T, inapplicable states are counted toward the sum of unique character-states.
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @param write If T, write the files locally.
#' @param log If T, write a file locally reporting the destination of each character.
#' @examples
#' # Synthetic example
#' splitNoStates(input = "../testdata/015_MORPH_data.nexus", input_format = "nexus", output_index = "../testdata/015_MORPH_data", ambiguity_addState = T, inapplicable_addState = T, log=T, write=T)
#' @references Černý, D., & Simonoff, A. L. (2023). Statistical evaluation of character support reveals the instability of higher-level dinosaur phylogeny. Scientific Reports, 13(1), 9273.
#' @references Khakurel, B., Grigsby, C., Tran, T. D., Zariwala, J., Höhna, S., & Wright, A. M. (2024). The fundamental role of character coding in Bayesian morphological phylogenetics. Systematic biology, 73(5), 861-871.
#' @export
splitNoStates <- function(input,
                          input_format = NULL,
                          ambiguity_addState = FALSE,
                          inapplicable_addState = FALSE,
                          output_index = "output",
                          log = FALSE,
                          write = TRUE) {

  # 1. Load matrix
  if (is.matrix(input) || is.data.frame(input)) {
    mat <- as.matrix(input)

  } else {
    if (is.null(input_format)) {
      stop("If 'input' is a file path, you must provide input_format = 'nexus' or 'tnt'")
    }
    if (input_format == "nexus") {
      mat <- TreeTools::ReadCharacters(input)
    } else if (input_format == "tnt") {
      mat <- TreeTools::ReadTntCharacters(input)
    } else {
      stop("input_format must be 'nexus' or 'tnt'")
    }
    mat <- as.matrix(mat)
  }

  # 2. Count states per column
  if (ambiguity_addState) {

    # Simple count: everything counts as a unique state,
    # including ambiguity tokens and possibly "-"
    uniq_counts <- apply(mat, 2, function(x) length(unique(x)))

  } else {

    # Collapse ambiguities — extract only digits and optionally "-"
    count_states <- function(column) {

      # Split polystates, e.g., "01" → "0","1"
      states <- unlist(strsplit(column, split = ""))

      # Keep digits always
      digit_states <- states[grepl("[0-9]", states)]

      # Optionally include "-"
      if (inapplicable_addState) {
        dash_states <- states[states == "-"]
      } else {
        dash_states <- character(0)
      }

      # Combine allowed states
      keep_states <- c(digit_states, dash_states)

      length(unique(keep_states))
    }

    uniq_counts <- apply(mat, 2, count_states)
  }

  # 3. Split matrix
  state_sizes <- sort(unique(uniq_counts))
  message("Generating ", length(state_sizes), " split matrices: ",
          paste(state_sizes, collapse = ", "))

  mat_split <- lapply(state_sizes, function(k) {
    mat[, uniq_counts == k, drop = FALSE]
  })

  names(mat_split) <- paste0("mat", state_sizes)

  # 4. Write matrices locally
  if (write) {
    if (is.null(input_format)) {
      stop("To write files, input_format must be 'nexus' or 'tnt'")
    }

    for (k in state_sizes) {
      # Ouput file name
      submat <- mat_split[[paste0("mat", k)]]
      outfile <- paste0(output_index, k,
                        if (input_format == "nexus") ".nex" else ".tnt")

      if (input_format == "nexus") {
        # Write a temporary file
        write.nexus.data(submat, file=outfile, format="standard", interleaved=F)
        # Remove the interleave section to avoid errors in IQTREE
        temp = gsub("INTERLEAVE=NO", "", readLines(outfile))
        # Write a permanent file
        writeLines(temp, outfile)
      } else if (input_format == "tnt") {
        # TNT uses [01] for ambiguities; convert any (01)-style tokens before writing
        submat[] <- gsub("\\(([^)]+)\\)", "[\\1]", submat, perl = TRUE)
        TreeTools::WriteTntCharacters(submat, file = outfile)
      }

      message("Wrote file: ", outfile)
    }
  }

  # 5. Optional log
  if (log) {
    log_file <- paste0(output_index, "_log.txt")
    con <- file(log_file, open = "wt")

    writeLines("=== Split Matrix Log ===", con)
    writeLines(paste("Original matrix dimensions:",
                     nrow(mat), "x", ncol(mat)), con)
    writeLines(paste("ambiguity_addState:", ambiguity_addState), con)
    writeLines(paste("inapplicable_addState:", inapplicable_addState), con)
    writeLines("", con)

    for (k in state_sizes) {
      idx <- which(uniq_counts == k)
      writeLines(paste0("Matrix mat", k, " (", length(idx), " columns):"), con)
      writeLines(paste("  Original column indices:",
                       paste(idx, collapse = ", ")), con)
      writeLines("", con)
    }

    close(con)
    message("Log file written to: ", log_file)
  }

  return(mat_split)
}
