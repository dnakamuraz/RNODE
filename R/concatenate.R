#' @title concatenate
#' @name concatenate
#' @description Concatenate two phylogenetic matrices (morphological or molecular) by columns.
#'   Inputs can be local NEXUS or TNT files, or matrices already loaded in R.
#'   Taxa present in only one matrix are filled with missing data ('?').
#'   The input format and data types (morphology, DNA, protein) are detected automatically
#'   if not explicitly defined in headers.
#'   Also reports and translates ordered (additive) characters safely to output matrix arrays.
#' @author Daniel YM Nakamura
#' @param input1 First input: a local file path (NEXUS or TNT) or a matrix/data.frame already loaded in R.
#' @param input2 Second input: a local file path (NEXUS or TNT) or a matrix/data.frame already loaded in R.
#' @param input1_format Format of input1 file: 'nexus' or 'tnt'. Detected automatically if NULL (default).
#' @param input2_format Format of input2 file: 'nexus' or 'tnt'. Detected automatically if NULL (default).
#' @param output_file Output file path (e.g. "testdata/059_TE_data.tnt"). If NULL, the result is returned as an R object only.
#' @param output_format Output format: 'nexus' or 'tnt'. If NULL, inferred from the extension of output_file ('.tnt' -> tnt, otherwise nexus). Ignored when output_file is NULL.
#' @return Invisibly returns the concatenated matrix.
#' @examples
#' # Concatenate two TNT files
#' concatenate(
#'   input1 = "testdata/059_MORPH_data.tnt",
#'   input2 = "testdata/059_MOL_data.tnt",
#'   output_file = "testdata/059_TE_data.tnt",
#'   output_format = "tnt"
#' )
#'
#' # Return as R object without writing
#' mat <- concatenate(
#'   input1 = "testdata/059_MORPH_data.tnt",
#'   input2 = "testdata/059_MOL_data.tnt"
#' )
#'
#' @export
concatenate <- function(input1,
                        input2,
                        input1_format = NULL,
                        input2_format = NULL,
                        output_file = NULL,
                        output_format = NULL) {
  # ---------------------------------------------------------------------------
  # Internal helpers
  # ---------------------------------------------------------------------------

  detect_format_from_ext <- function(path) {
    if (grepl("\\.tnt$", path, ignore.case = TRUE)) "tnt" else "nexus"
  }

  detect_type_from_nexus <- function(nexus_file) {
    lines <- tolower(readLines(nexus_file))
    fmt <- grep("\\bformat\\b", lines, value = TRUE)
    if (length(fmt) > 0) {
      if (grepl("datatype\\s*=\\s*protein|\\bprotein\\b", fmt[1])) {
        return("protein")
      }
      if (grepl("datatype\\s*=\\s*dna|\\bdna\\b", fmt[1])) {
        return("dna")
      }
      if (grepl("standard", fmt[1])) {
        return("standard")
      }
      if (grepl("mixed", fmt[1])) {
        return("mixed")
      }
    }
    return("unknown") # Fallback to matrix character detection if no format found
  }

  detect_type_from_tnt <- function(tnt_file) {
    lines <- readLines(tnt_file)
    lines_lc <- tolower(lines)
    if (any(grepl("nstates\\s+(prot|32)", lines_lc))) {
      return("protein")
    }
    if (any(grepl("&\\s*\\[prot\\]", lines_lc))) {
      return("protein")
    }
    if (any(grepl("nstates\\s+dna", lines_lc))) {
      return("dna")
    }
    if (any(grepl("&\\s*\\[dna\\]", lines_lc))) {
      return("dna")
    }
    return("unknown") # Fallback to matrix character detection if headless
  }

  # Helper detecting pure arrays from unflagged inputs via alphabet signatures
  detect_matrix_type <- function(mat) {
    vals <- unique(toupper(as.character(mat)))
    vals <- setdiff(vals, c("?", "-", "N", "X", "(", ")", "[", "]", "{", "}"))
    if (length(vals) == 0) {
      return("standard")
    }

    # If there are numerical states (0-9), it is inherently morphology/standard
    if (any(grepl("[0-9]", vals))) {
      return("standard")
    }

    # Check for protein-exclusive identifiers first to prevent misfires
    # (letters absent from normal ambiguous nucleotide dictionaries)
    aa_exclusive <- c("E", "F", "I", "L", "P", "Q", "Z", "J")
    if (any(vals %in% aa_exclusive)) {
      return("protein")
    }

    # Strict IUPAC/ambiguous DNA bounds
    dna_chars <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")
    if (all(vals %in% dna_chars)) {
      return("dna")
    }

    # Strict IUPAC protein bounds
    aa_chars <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J")
    if (all(vals %in% aa_chars)) {
      return("protein")
    }

    # Mixed/ambiguous characters defaults back to morphology array
    return("standard")
  }

  get_ordered_chars <- function(input_val, format_type) {
    if (is.matrix(input_val) || is.data.frame(input_val) || !is.character(input_val) || !file.exists(input_val)) {
      return(integer(0))
    }
    lines <- readLines(input_val, warn = FALSE)
    content <- paste(lines, collapse = " ")
    ordered_idx <- integer(0)

    if (format_type == "nexus") {
      matches <- gregexpr("(?i)TYPESET\\s+\\*?\\s*Ordered\\s*=\\s*([^;]+);", content, perl = TRUE)
      if (matches[[1]][1] != -1) {
        captures <- regmatches(content, matches)[[1]]
        for (cap in captures) {
          body <- sub("(?i).*?Ordered\\s*=\\s*", "", cap, perl = TRUE)
          body <- sub(";.*$", "", body)
          tokens <- unlist(strsplit(trimws(body), "\\s+"))
          for (tk in tokens) {
            if (grepl("-", tk)) {
              pts <- as.numeric(unlist(strsplit(tk, "-")))
              if (length(pts) == 2 && !any(is.na(pts))) ordered_idx <- c(ordered_idx, seq(pts[1], pts[2]))
            } else if (tk != "") {
              v <- as.numeric(tk)
              if (!is.na(v)) ordered_idx <- c(ordered_idx, v)
            }
          }
        }
      }
      return(unique(ordered_idx))
    } else if (format_type == "tnt") {
      matches <- gregexpr("(?i)cc\\s+\\+\\s*([^;]+);", content, perl = TRUE)
      if (matches[[1]][1] != -1) {
        captures <- regmatches(content, matches)[[1]]
        for (cap in captures) {
          body <- sub("(?i)cc\\s+\\+\\s*", "", cap, perl = TRUE)
          body <- sub(";.*$", "", body)
          tokens <- unlist(strsplit(trimws(body), "\\s+"))
          for (tk in tokens) {
            if (grepl("\\.", tk)) {
              parts <- as.numeric(strsplit(tk, "\\.+")[[1]])
              if (length(parts) == 2 && !any(is.na(parts))) {
                ordered_idx <- c(ordered_idx, seq(parts[1] + 1, parts[2] + 1))
              }
            } else if (grepl("-", tk)) {
              parts <- as.numeric(strsplit(tk, "-")[[1]])
              if (length(parts) == 2 && !any(is.na(parts))) {
                ordered_idx <- c(ordered_idx, seq(parts[1] + 1, parts[2] + 1))
              }
            } else if (tk != "") {
              v <- as.numeric(tk)
              if (!is.na(v)) ordered_idx <- c(ordered_idx, v + 1)
            }
          }
        }
      }
      return(unique(ordered_idx))
    }
    return(integer(0))
  }

  format_tnt_ordered <- function(idx) {
    if (length(idx) == 0) {
      return("")
    }
    idx <- sort(unique(idx))
    chunks <- character(0)
    start_val <- prev_val <- idx[1]
    if (length(idx) > 1) {
      for (i in seq_along(idx)[-1]) {
        if (idx[i] == prev_val + 1) {
          prev_val <- idx[i]
        } else {
          chunks <- c(chunks, if (start_val == prev_val) start_val else paste0(start_val, ".", prev_val))
          start_val <- prev_val <- idx[i]
        }
      }
    }
    chunks <- c(chunks, if (start_val == prev_val) start_val else paste0(start_val, ".", prev_val))
    paste(chunks, collapse = " ")
  }

  format_nexus_ordered <- function(idx) {
    if (length(idx) == 0) {
      return("")
    }
    idx <- sort(unique(idx))
    chunks <- character(0)
    start_val <- prev_val <- idx[1]
    if (length(idx) > 1) {
      for (i in seq_along(idx)[-1]) {
        if (idx[i] == prev_val + 1) {
          prev_val <- idx[i]
        } else {
          chunks <- c(chunks, if (start_val == prev_val) start_val else paste0(start_val, "-", prev_val))
          start_val <- prev_val <- idx[i]
        }
      }
    }
    chunks <- c(chunks, if (start_val == prev_val) start_val else paste0(start_val, "-", prev_val))
    paste(chunks, collapse = " ")
  }

  read_tnt_manual <- function(tnt_file) {
    lines <- readLines(tnt_file)
    xread_idx <- which(grepl("^xread", lines, ignore.case = TRUE))[1]
    if (is.na(xread_idx)) stop("Could not find 'xread' block in TNT file: ", tnt_file)
    lines <- lines[xread_idx:length(lines)]
    header_parts <- strsplit(trimws(lines[2]), "\\s+")[[1]]
    nchars <- as.numeric(header_parts[1])
    ntaxa <- as.numeric(header_parts[2])
    if (is.na(nchars) || is.na(ntaxa)) stop("Could not parse TNT dimensions in file: ", tnt_file)

    data_lines <- trimws(lines[-c(1, 2)])
    data_lines <- data_lines[nzchar(data_lines)]
    data_lines <- data_lines[!grepl("^&\\s*\\[", data_lines)]
    end_idx <- which(data_lines == ";")[1]
    if (!is.na(end_idx)) data_lines <- data_lines[seq_len(end_idx - 1)]
    if (length(data_lines) < ntaxa) stop("Expected ", ntaxa, " taxa but found ", length(data_lines), " sequence lines in: ", tnt_file)

    seq_lines <- data_lines[seq_len(ntaxa)]
    seq_parts <- regmatches(seq_lines, regexec("^([^\\s]+)\\s+(.+)$", seq_lines, perl = TRUE))
    parsed <- vapply(seq_parts, length, integer(1)) == 3
    if (!all(parsed)) stop("Could not parse taxon/sequence on line(s): ", paste(which(!parsed), collapse = ", "))

    taxa_names <- vapply(seq_parts, `[`, character(1), 2)
    sequences <- vapply(seq_parts, `[`, character(1), 3)
    sequences <- gsub("\\s+", "", sequences, perl = TRUE)

    if (!any(grepl("[\\[\\{\\(]", sequences, perl = TRUE))) {
      seq_widths <- nchar(sequences, type = "chars")
      parsed_nchars <- max(seq_widths)
      if (!all(seq_widths == parsed_nchars)) {
        split_sequences <- strsplit(sequences, "", fixed = TRUE)
        split_sequences <- lapply(split_sequences, function(x) {
          if (length(x) < parsed_nchars) c(x, rep("?", parsed_nchars - length(x))) else x
        })
        mat <- matrix(unlist(split_sequences, use.names = FALSE), nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
        rownames(mat) <- taxa_names
        warning("TNT sequences have unequal lengths; shorter sequences were padded with '?'.")
        return(mat)
      }
      mat <- matrix(unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE), nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
      rownames(mat) <- taxa_names
      return(mat)
    }

    token_pattern <- "\\[[^]]+\\]|\\{[^}]+\\}|\\([^)]*\\)|."
    mat_list <- regmatches(sequences, gregexpr(token_pattern, sequences, perl = TRUE))
    token_counts <- lengths(mat_list)
    parsed_nchars <- max(token_counts)
    if (!all(token_counts == parsed_nchars)) {
      mat_list <- lapply(mat_list, function(x) {
        if (length(x) < parsed_nchars) c(x, rep("?", parsed_nchars - length(x))) else x
      })
      warning("TNT sequences have unequal token counts; shorter sequences were padded with '?'.")
    }
    mat <- matrix(unlist(mat_list, use.names = FALSE), nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
    rownames(mat) <- taxa_names
    return(mat)
  }

  tnt_fix_ambiguity <- function(sequence) gsub("\\(([^)]+)\\)", "[\\1]", sequence, perl = TRUE)

  write_tnt_lines <- function(lines, filename) {
    if (length(lines) > 0 && lines[1] == "nstates 32") {
      text <- paste0(lines[1], "\n", paste(lines[-1], collapse = "\r\n"), "\r\n")
      writeBin(charToRaw(text), filename)
    } else {
      writeLines(lines, filename, sep = "\r\n")
    }
  }

  write_tnt_concat <- function(mat1, mat2, filename, type1 = "standard", type2 = "standard", ordered1 = integer(0), ordered2 = integer(0)) {
    ntaxa <- nrow(mat1)
    total_nchars <- ncol(mat1) + ncol(mat2)

    lines <- c("xread", paste(total_nchars, ntaxa))

    # Identify matching block headers explicitly for appending sequentially
    header_of <- function(t) {
      switch(t,
        protein = "&[prot]",
        dna = "&[dna]",
        "&[num]"
      )
    }

    write_block <- function(m, h) {
      lines <<- c(lines, "", h)
      for (i in seq_len(ntaxa)) {
        seq <- tnt_fix_ambiguity(paste(m[i, ], collapse = ""))
        lines <<- c(lines, paste0(rownames(m)[i], "\t", seq))
      }
    }

    write_block(mat1, header_of(type1))
    write_block(mat2, header_of(type2))

    # Safely merge zero-based tracking adjustments for explicit TNT arrays
    combined_tnt_ordered <- integer(0)
    if (length(ordered1) > 0) combined_tnt_ordered <- c(combined_tnt_ordered, ordered1 - 1)
    if (length(ordered2) > 0) combined_tnt_ordered <- c(combined_tnt_ordered, ordered2 - 1 + ncol(mat1))

    lines <- c(lines, ";", "")
    if (length(combined_tnt_ordered) > 0) {
      lines <- c(lines, paste0("cc + ", format_tnt_ordered(combined_tnt_ordered), " ;"), "")
    }

    lines <- c(lines, "proc /;", "comments 0", ";", "", "")
    write_tnt_lines(lines, filename)
  }

  # ---------------------------------------------------------------------------
  # Load input1
  # ---------------------------------------------------------------------------
  ordered1 <- integer(0)
  type1 <- "standard"
  if (is.matrix(input1) || is.data.frame(input1)) {
    mat1 <- as.matrix(input1)
    type1 <- detect_matrix_type(mat1)
  } else {
    if (is.null(input1_format)) input1_format <- detect_format_from_ext(input1)
    input1_format <- match.arg(tolower(input1_format), c("nexus", "tnt"))
    ordered1 <- get_ordered_chars(input1, input1_format)

    if (input1_format == "nexus") {
      header_type1 <- detect_type_from_nexus(input1)
      mat1 <- as.matrix(TreeTools::ReadCharacters(input1))
    } else {
      header_type1 <- detect_type_from_tnt(input1)
      mat1 <- read_tnt_manual(input1)
    }

    type1 <- if (header_type1 != "unknown") header_type1 else detect_matrix_type(mat1)
  }

  # ---------------------------------------------------------------------------
  # Load input2
  # ---------------------------------------------------------------------------
  ordered2 <- integer(0)
  type2 <- "standard"
  if (is.matrix(input2) || is.data.frame(input2)) {
    mat2 <- as.matrix(input2)
    type2 <- detect_matrix_type(mat2)
  } else {
    if (is.null(input2_format)) input2_format <- detect_format_from_ext(input2)
    input2_format <- match.arg(tolower(input2_format), c("nexus", "tnt"))
    ordered2 <- get_ordered_chars(input2, input2_format)

    if (input2_format == "nexus") {
      header_type2 <- detect_type_from_nexus(input2)
      mat2 <- as.matrix(TreeTools::ReadCharacters(input2))
    } else {
      header_type2 <- detect_type_from_tnt(input2)
      mat2 <- read_tnt_manual(input2)
    }

    type2 <- if (header_type2 != "unknown") header_type2 else detect_matrix_type(mat2)
  }

  # ---------------------------------------------------------------------------
  # Validate & Merge Matrices
  # ---------------------------------------------------------------------------
  taxa1 <- rownames(mat1)
  taxa2 <- rownames(mat2)
  if (is.null(taxa1) || is.null(taxa2)) stop("Both matrices must have row names (taxon names).")

  all_taxa <- union(taxa1, taxa2)

  cat("Taxa in matrix 1       :", length(taxa1), "\n")
  cat("Taxa in matrix 2       :", length(taxa2), "\n")
  cat("Shared taxa            :", length(intersect(taxa1, taxa2)), "\n")
  cat("Total taxa in output   :", length(all_taxa), "\n")

  cat("Characters in matrix 1 :", ncol(mat1))
  if (length(ordered1) > 0) cat(sprintf(" (including %d ordered/additive)", length(ordered1)))
  cat(sprintf(" - Detected Block [%s]", toupper(type1)))
  cat("\n")

  cat("Characters in matrix 2 :", ncol(mat2))
  if (length(ordered2) > 0) cat(sprintf(" (including %d ordered/additive)", length(ordered2)))
  cat(sprintf(" - Detected Block [%s]", toupper(type2)))
  cat("\n")

  # Standard shift (Appending rightwise against total nexus base)
  combined_ordered_nexus <- unique(c(
    if (length(ordered1) > 0) ordered1 else integer(0),
    if (length(ordered2) > 0) ordered2 + ncol(mat1) else integer(0)
  ))

  cat("Total characters       :", ncol(mat1) + ncol(mat2))
  if (length(combined_ordered_nexus) > 0) cat(sprintf(" (including %d ordered/additive)", length(combined_ordered_nexus)))
  cat("\n\n")

  align_to_taxa <- function(mat, taxa) {
    out <- matrix("?", nrow = length(taxa), ncol = ncol(mat), dimnames = list(taxa, colnames(mat)))
    present <- intersect(taxa, rownames(mat))
    out[present, ] <- mat[present, , drop = FALSE]
    out
  }

  mat1_aligned <- align_to_taxa(mat1, all_taxa)
  mat2_aligned <- align_to_taxa(mat2, all_taxa)

  # Column-Wise Nexus append internally mapped natively
  result <- cbind(mat1_aligned, mat2_aligned)
  if (length(combined_ordered_nexus) > 0) attr(result, "ordered") <- combined_ordered_nexus
  attr(result, "type1") <- type1
  attr(result, "type2") <- type2

  # ---------------------------------------------------------------------------
  # Write output file
  # ---------------------------------------------------------------------------
  if (!is.null(output_file)) {
    if (is.null(output_format)) {
      output_format <- if (grepl("\\.tnt$", output_file, ignore.case = TRUE)) "tnt" else "nexus"
    }
    output_format <- match.arg(tolower(output_format), c("nexus", "tnt"))

    out_dir <- dirname(output_file)
    if (out_dir != "." && !dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir, call. = FALSE)

    if (output_format == "tnt") {
      if (!grepl("\\.tnt$", output_file, ignore.case = TRUE)) output_file <- paste0(output_file, ".tnt")

      write_tnt_concat(mat1_aligned, mat2_aligned, output_file, type1, type2, ordered1, ordered2)

      cat("Concatenated matrix written to:", output_file, "\n")
      cat("Format: TNT (blocks:", type1, "+", type2, ")\n")
    } else {
      if (!grepl("\\.(nexus|nex)$", output_file, ignore.case = TRUE)) output_file <- paste0(output_file, ".nexus")

      # Using format 'standard' dynamically preserves any matrix layout gracefully across variables (numbers, A, C)
      ape::write.nexus.data(result, file = output_file, format = "standard", interleaved = FALSE)

      tmp <- readLines(output_file)
      tmp <- gsub("INTERLEAVE=NO", "", tmp)
      tmp <- gsub("write.nexus.data.R", "RNODE", tmp)

      assump_block <- character(0)

      # Append DATATYPE format explicitly translating Mixed Arrays & Layout coordinates to standards blocks
      if (type1 != type2) {
        mixed_str <- paste0("DATATYPE=MIXED(", toupper(type1), ":1-", ncol(mat1), ", ", toupper(type2), ":", ncol(mat1) + 1, "-", ncol(mat1) + ncol(mat2), ")")
        tmp <- sub("DATATYPE=\\w+", mixed_str, tmp, ignore.case = TRUE)

        # Output Interleaved Nexus Sets Header defining standard partition arrays
        assump_block <- c(
          assump_block,
          "BEGIN SETS;",
          paste0("    CHARSET ", type1, "_block = 1-", ncol(mat1), ";"),
          paste0("    CHARSET ", type2, "_block = ", ncol(mat1) + 1, "-", ncol(mat1) + ncol(mat2), ";"),
          "END;"
        )
      } else {
        tmp <- sub("DATATYPE=\\w+", paste0("DATATYPE=", toupper(type1)), tmp, ignore.case = TRUE)
      }

      if (length(combined_ordered_nexus) > 0) {
        assump_block <- c(
          assump_block,
          "BEGIN ASSUMPTIONS;",
          "    OPTIONS DEFTYPE=unord;",
          paste0("    TYPESET * Ordered = ", format_nexus_ordered(combined_ordered_nexus), ";"),
          "END;"
        )
      }

      if (length(assump_block) > 0) tmp <- c(tmp, assump_block)

      writeLines(tmp, output_file)
      cat("Concatenated matrix written to:", output_file, "\n")
      cat("Format: NEXUS", if (type1 != type2) paste0("(MIXED ", type1, "+", type2, ")") else paste0("(", type1, ")"), "\n")
    }
  }

  invisible(result)
}
