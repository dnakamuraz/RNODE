#' @title concatenate
#' @name concatenate
#' @description Concatenate two phylogenetic matrices (morphological or molecular) by columns.
#'   Inputs can be local NEXUS or TNT files, or matrices already loaded in R.
#'   Taxa present in only one matrix are filled with missing data ('?').
#'   The input format is detected automatically when not provided.
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
#' concatenate(input1 = "testdata/059_MORPH_data.tnt",
#'             input2 = "testdata/059_MOL_data.tnt",
#'             output_file = "testdata/059_TE_data.tnt",
#'             output_format = "tnt")
#'
#' # Return as R object without writing
#' mat <- concatenate(input1 = "testdata/059_MORPH_data.tnt",
#'                    input2 = "testdata/059_MOL_data.tnt")
#'
#' @export
concatenate <- function(input1,
                        input2,
                        input1_format = NULL,
                        input2_format = NULL,
                        output_file   = NULL,
                        output_format = NULL) {

  # ---------------------------------------------------------------------------
  # Internal helpers
  # ---------------------------------------------------------------------------

  # Auto-detect format from file extension
  detect_format_from_ext <- function(path) {
    if (grepl("\\.tnt$", path, ignore.case = TRUE)) "tnt" else "nexus"
  }

  # Detect data type (standard / dna / protein) from a NEXUS file
  detect_type_from_nexus <- function(nexus_file) {
    lines <- tolower(readLines(nexus_file))
    fmt   <- grep("\\bformat\\b", lines, value = TRUE)
    if (length(fmt) > 0) {
      if (grepl("datatype\\s*=\\s*protein|protein", fmt[1])) return("protein")
      if (grepl("datatype\\s*=\\s*dna|dna",         fmt[1])) return("dna")
      if (grepl("standard",                          fmt[1])) return("standard")
    }
    "standard"
  }

  # Detect data type from a TNT file
  detect_type_from_tnt <- function(tnt_file) {
    lines <- tolower(readLines(tnt_file, n = 10))
    if (any(grepl("nstates\\s+(prot|32)", lines))) return("protein")
    "standard"
  }

  # Read a TNT file, returning a character matrix (taxa x characters)
  read_tnt_manual <- function(tnt_file) {
    lines <- readLines(tnt_file)

    xread_idx <- which(grepl("^xread", lines, ignore.case = TRUE))[1]
    if (is.na(xread_idx)) {
      stop("Could not find 'xread' block in TNT file: ", tnt_file)
    }
    lines <- lines[xread_idx:length(lines)]

    header_parts <- strsplit(trimws(lines[2]), "\\s+")[[1]]
    nchars <- as.numeric(header_parts[1])
    ntaxa  <- as.numeric(header_parts[2])
    if (is.na(nchars) || is.na(ntaxa)) {
      stop("Could not parse TNT dimensions in file: ", tnt_file)
    }

    data_lines <- trimws(lines[-c(1, 2)])
    data_lines <- data_lines[nzchar(data_lines)]
    data_lines <- data_lines[!grepl("^&\\s*\\[", data_lines)]
    end_idx    <- which(data_lines == ";")[1]
    if (!is.na(end_idx)) data_lines <- data_lines[seq_len(end_idx - 1)]

    if (length(data_lines) < ntaxa) {
      stop("Expected ", ntaxa, " taxa but found ", length(data_lines),
           " sequence lines in: ", tnt_file)
    }
    seq_lines <- data_lines[seq_len(ntaxa)]

    seq_parts <- regmatches(
      seq_lines,
      regexec("^([^\\s]+)\\s+(.+)$", seq_lines, perl = TRUE)
    )
    parsed <- vapply(seq_parts, length, integer(1)) == 3
    if (!all(parsed)) {
      stop("Could not parse taxon/sequence on line(s): ",
           paste(which(!parsed), collapse = ", "))
    }

    taxa_names <- vapply(seq_parts, `[`, character(1), 2)
    sequences  <- vapply(seq_parts, `[`, character(1), 3)
    sequences  <- gsub("\\s+", "", sequences, perl = TRUE)

    # Fast path: no polymorphic tokens
    if (!any(grepl("[\\[\\{\\(]", sequences, perl = TRUE))) {
      seq_widths    <- nchar(sequences, type = "chars")
      parsed_nchars <- max(seq_widths)
      if (!all(seq_widths == parsed_nchars)) {
        split_sequences <- strsplit(sequences, "", fixed = TRUE)
        split_sequences <- lapply(split_sequences, function(x) {
          if (length(x) < parsed_nchars) c(x, rep("?", parsed_nchars - length(x))) else x
        })
        mat <- matrix(unlist(split_sequences, use.names = FALSE),
                      nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
        rownames(mat) <- taxa_names
        warning("TNT sequences have unequal lengths; shorter sequences were padded with '?'.")
        return(mat)
      }
      mat <- matrix(
        unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE),
        nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE
      )
      rownames(mat) <- taxa_names
      return(mat)
    }

    # Slow path: polymorphic / inapplicable tokens
    token_pattern <- "\\[[^]]+\\]|\\{[^}]+\\}|\\([^)]*\\)|."
    mat_list      <- regmatches(sequences, gregexpr(token_pattern, sequences, perl = TRUE))
    token_counts  <- lengths(mat_list)
    parsed_nchars <- max(token_counts)
    if (!all(token_counts == parsed_nchars)) {
      mat_list <- lapply(mat_list, function(x) {
        if (length(x) < parsed_nchars) c(x, rep("?", parsed_nchars - length(x))) else x
      })
      warning("TNT sequences have unequal token counts; shorter sequences were padded with '?'.")
    }
    mat <- matrix(unlist(mat_list, use.names = FALSE),
                  nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
    rownames(mat) <- taxa_names
    return(mat)
  }

  # Convert (01)-style ambiguity tokens to TNT [01] notation in a sequence string
  tnt_fix_ambiguity <- function(sequence) {
    gsub("\\(([^)]+)\\)", "[\\1]", sequence, perl = TRUE)
  }

  # Write lines to a TNT file (CRLF endings; nstates 32 header gets LF separator)
  write_tnt_lines <- function(lines, filename) {
    if (length(lines) > 0 && lines[1] == "nstates 32") {
      text <- paste0(lines[1], "\n", paste(lines[-1], collapse = "\r\n"), "\r\n")
      writeBin(charToRaw(text), filename)
    } else {
      writeLines(lines, filename, sep = "\r\n")
    }
  }

  # Build and write a single-block TNT file from a matrix
  write_tnt_single <- function(mat, filename, data_type = "standard") {
    ntaxa  <- nrow(mat)
    nchars <- ncol(mat)

    lines <- c("xread", paste(nchars, ntaxa))
    if (data_type == "protein") lines <- c("nstates 32", lines)

    block_header <- switch(data_type,
                           protein  = "& [prot]",
                           dna      = "& [dna]",
                           "& [num]")
    lines <- c(lines, "", block_header)

    for (i in seq_len(ntaxa)) {
      seq <- tnt_fix_ambiguity(paste(mat[i, ], collapse = ""))
      lines <- c(lines, paste0(rownames(mat)[i], "\t", seq))
    }
    lines <- c(lines, ";", "", "", "proc /;", "comments 0", ";", "", "")
    write_tnt_lines(lines, filename)
  }

  # Build and write a two-block TNT file from two matrices (concatenated)
  write_tnt_concat <- function(mat1, mat2, filename, type1 = "standard", type2 = "standard") {
    ntaxa        <- nrow(mat1)
    total_nchars <- ncol(mat1) + ncol(mat2)

    lines <- c("xread", paste(total_nchars, ntaxa))
    if (type1 == "protein" || type2 == "protein") lines <- c("nstates 32", lines)

    header_of <- function(t) switch(t, protein = "& [prot]", dna = "& [dna]", "& [num]")

    # Protein block always comes first in TNT convention
    write_block <- function(m, h) {
      lines <<- c(lines, "", h)
      for (i in seq_len(ntaxa)) {
        seq   <- tnt_fix_ambiguity(paste(m[i, ], collapse = ""))
        lines <<- c(lines, paste0(rownames(m)[i], "\t", seq))
      }
    }

    if (type1 == "protein" && type2 != "protein") {
      write_block(mat1, header_of(type1))
      write_block(mat2, header_of(type2))
    } else if (type1 != "protein" && type2 == "protein") {
      write_block(mat2, header_of(type2))
      write_block(mat1, header_of(type1))
    } else {
      write_block(mat1, header_of(type1))
      write_block(mat2, header_of(type2))
    }

    lines <- c(lines, ";", "", "", "proc /;", "comments 0", ";", "", "")
    write_tnt_lines(lines, filename)
  }

  # ---------------------------------------------------------------------------
  # Load input1
  # ---------------------------------------------------------------------------
  type1 <- "standard"

  if (is.matrix(input1) || is.data.frame(input1)) {
    mat1 <- as.matrix(input1)
  } else {
    if (is.null(input1_format)) input1_format <- detect_format_from_ext(input1)
    input1_format <- match.arg(tolower(input1_format), c("nexus", "tnt"))

    if (input1_format == "nexus") {
      type1 <- detect_type_from_nexus(input1)
      mat1  <- as.matrix(TreeTools::ReadCharacters(input1))
    } else {
      type1 <- detect_type_from_tnt(input1)
      mat1  <- read_tnt_manual(input1)
    }
  }

  # ---------------------------------------------------------------------------
  # Load input2
  # ---------------------------------------------------------------------------
  type2 <- "standard"

  if (is.matrix(input2) || is.data.frame(input2)) {
    mat2 <- as.matrix(input2)
  } else {
    if (is.null(input2_format)) input2_format <- detect_format_from_ext(input2)
    input2_format <- match.arg(tolower(input2_format), c("nexus", "tnt"))

    if (input2_format == "nexus") {
      type2 <- detect_type_from_nexus(input2)
      mat2  <- as.matrix(TreeTools::ReadCharacters(input2))
    } else {
      type2 <- detect_type_from_tnt(input2)
      mat2  <- read_tnt_manual(input2)
    }
  }

  # ---------------------------------------------------------------------------
  # Validate row names
  # ---------------------------------------------------------------------------
  taxa1 <- rownames(mat1)
  taxa2 <- rownames(mat2)
  if (is.null(taxa1) || is.null(taxa2)) {
    stop("Both matrices must have row names (taxon names).")
  }

  # Union of all taxa, preserving order (mat1 first, then mat2-unique)
  all_taxa <- union(taxa1, taxa2)

  cat("Taxa in matrix 1       :", length(taxa1), "\n")
  cat("Taxa in matrix 2       :", length(taxa2), "\n")
  cat("Shared taxa            :", length(intersect(taxa1, taxa2)), "\n")
  cat("Total taxa in output   :", length(all_taxa), "\n")
  cat("Characters in matrix 1 :", ncol(mat1), "\n")
  cat("Characters in matrix 2 :", ncol(mat2), "\n")
  cat("Total characters       :", ncol(mat1) + ncol(mat2), "\n\n")

  # ---------------------------------------------------------------------------
  # Align both matrices to the full taxon set (fill absent rows with '?')
  # ---------------------------------------------------------------------------
  align_to_taxa <- function(mat, taxa) {
    out <- matrix("?", nrow = length(taxa), ncol = ncol(mat),
                  dimnames = list(taxa, colnames(mat)))
    present <- intersect(taxa, rownames(mat))
    out[present, ] <- mat[present, , drop = FALSE]
    out
  }

  mat1_aligned <- align_to_taxa(mat1, all_taxa)
  mat2_aligned <- align_to_taxa(mat2, all_taxa)

  # Concatenate column-wise
  result <- cbind(mat1_aligned, mat2_aligned)

  # ---------------------------------------------------------------------------
  # Write output file (if requested)
  # ---------------------------------------------------------------------------
  if (!is.null(output_file)) {

    # Determine output format
    if (is.null(output_format)) {
      output_format <- if (grepl("\\.tnt$", output_file, ignore.case = TRUE)) "tnt" else "nexus"
    }
    output_format <- match.arg(tolower(output_format), c("nexus", "tnt"))

    # Validate output directory
    out_dir <- dirname(output_file)
    if (out_dir != "." && !dir.exists(out_dir)) {
      stop("Output directory does not exist: ", out_dir, call. = FALSE)
    }

    if (output_format == "tnt") {
      # Ensure .tnt extension
      if (!grepl("\\.tnt$", output_file, ignore.case = TRUE)) {
        output_file <- paste0(output_file, ".tnt")
      }
      write_tnt_concat(mat1_aligned, mat2_aligned, output_file, type1, type2)
      cat("Concatenated matrix written to:", output_file, "\n")
      cat("Format: TNT (blocks:", type1, "+", type2, ")\n")

    } else {
      # Ensure .nexus extension
      if (!grepl("\\.(nexus|nex)$", output_file, ignore.case = TRUE)) {
        output_file <- paste0(output_file, ".nexus")
      }
      nexus_type <- if (type1 == "protein" || type2 == "protein") {
        "protein"
      } else if (type1 == "dna" || type2 == "dna") {
        "dna"
      } else {
        "standard"
      }
      ape::write.nexus.data(result, file = output_file,
                            format = nexus_type, interleaved = FALSE)
      tmp <- readLines(output_file)
      tmp <- gsub("INTERLEAVE=NO", "", tmp)
      tmp <- gsub("write.nexus.data.R", "RNODE", tmp)
      writeLines(tmp, output_file)
      cat("Concatenated matrix written to:", output_file, "\n")
      cat("Format: NEXUS (", nexus_type, ")\n")
    }
  }

  invisible(result)
}
