#' @title concatenate
#' @name concatenate
#' @description Concatenate two or more phylogenetic matrices (morphological or molecular) by columns.
#' Inputs can be local NEXUS or TNT files, or matrices already loaded in R.
#' Taxa present in only one matrix are filled with missing data ('?').
#' The input format and data types (morphology, DNA, protein) are detected automatically
#' if not explicitly defined in headers.
#' Also reports and translates ordered (additive) characters safely to output matrix arrays.
#' @author Daniel YM Nakamura
#' @param input1 First input: a local file path (NEXUS or TNT) or a matrix/data.frame already loaded in R. Ignored if `multiple_input` is provided.
#' @param input2 Second input: a local file path (NEXUS or TNT) or a matrix/data.frame already loaded in R. Ignored if `multiple_input` is provided.
#' @param input1_format Format of input1 file: 'nexus' or 'tnt'. Detected automatically if NULL (default).
#' @param input2_format Format of input2 file: 'nexus' or 'tnt'. Detected automatically if NULL (default).
#' @param output_file Output file path (e.g. "testdata/059_TE_data.tnt"). If NULL, the result is returned as an R object only.
#' @param output_format Output format: 'nexus' or 'tnt'. If NULL, inferred from the extension of output_file ('.tnt' -> tnt, otherwise nexus). Ignored when output_file is NULL.
#' @param multiple_input A directory path containing multiple TNT/NEXUS files to concatenate. If provided, `input1` and `input2` are ignored.
#' @return Invisibly returns the concatenated matrix.
#' @examples
#' # Concatenate all files in a folder
#' concatenate(multiple_input = "testdata/partitioned/", output_file = "combined.tnt")
#' @export
concatenate <- function(input1 = NULL, input2 = NULL,
                        input1_format = NULL, input2_format = NULL,
                        output_file = NULL, output_format = NULL,
                        multiple_input = NULL) {
  # ---------------------------------------------------------------------------
  # Internal helpers (unchanged from previous version)
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
    return("unknown")
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
    return("unknown")
  }

  detect_matrix_type <- function(mat) {
    vals <- unique(toupper(as.character(mat)))
    vals <- setdiff(vals, c("?", "-", "N", "X", "(", ")", "[", "]", "{", "}"))
    if (length(vals) == 0) {
      return("standard")
    }
    if (any(grepl("[0-9]", vals))) {
      return("standard")
    }
    aa_exclusive <- c("E", "F", "I", "L", "P", "Q", "Z", "J")
    if (any(vals %in% aa_exclusive)) {
      return("protein")
    }
    dna_chars <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")
    if (all(vals %in% dna_chars)) {
      return("dna")
    }
    aa_chars <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J")
    if (all(vals %in% aa_chars)) {
      return("protein")
    }
    return("standard")
  }

  get_ordered_chars <- function(input_val, format_type) {
    if (is.matrix(input_val) || is.data.frame(input_val) ||
      !is.character(input_val) || !file.exists(input_val)) {
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
              if (length(parts) == 2 && !any(is.na(parts))) ordered_idx <- c(ordered_idx, seq(parts[1] + 1, parts[2] + 1))
            } else if (grepl("-", tk)) {
              parts <- as.numeric(strsplit(tk, "-")[[1]])
              if (length(parts) == 2 && !any(is.na(parts))) ordered_idx <- c(ordered_idx, seq(parts[1] + 1, parts[2] + 1))
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

  format_tnt_sequence <- function(cells) {
    sapply(cells, function(x) {
      if (grepl("^[\\[\\(\\{].*[\\]\\)\\}]$", x)) {
        return(x)
      }
      if (nchar(x) > 1) {
        return(paste0("[", x, "]"))
      }
      return(x)
    }, USE.NAMES = FALSE)
  }

  # Robust TNT reader (handles titles, interleaved blocks, brackets)
  read_tnt_manual <- function(tnt_file) {
    lines <- readLines(tnt_file)
    xread_idx <- which(grepl("^xread", lines, ignore.case = TRUE))[1]
    if (is.na(xread_idx)) stop("Could not find 'xread' block in TNT file: ", tnt_file)

    extract_dims <- function(txt) {
      txt <- gsub("'[^']*'", "", txt)
      nums <- regmatches(txt, gregexpr("\\d+", txt))[[1]]
      if (length(nums) >= 2) {
        return(list(nchars = as.numeric(nums[1]), ntaxa = as.numeric(nums[2])))
      }
      return(NULL)
    }

    dims <- extract_dims(lines[xread_idx])
    data_start <- xread_idx + 1
    if (is.null(dims)) {
      for (k in (xread_idx + 1):min(xread_idx + 3, length(lines))) {
        dims <- extract_dims(lines[k])
        if (!is.null(dims)) {
          data_start <- k + 1
          break
        }
      }
      if (is.null(dims)) stop("Cannot find dimensions (nchars ntaxa) after xread in: ", tnt_file)
    }
    nchars <- dims$nchars
    ntaxa <- dims$ntaxa

    data_lines <- lines[data_start:length(lines)]
    end_idx <- which(grepl("^\\s*;", data_lines))[1]
    if (!is.na(end_idx)) data_lines <- data_lines[1:(end_idx - 1)]

    taxa_order <- character(ntaxa)
    seq_accum <- character(ntaxa)
    block_sizes <- integer(0)
    block_types <- character(0)
    block_taxon_idx <- 1
    current_block_seq <- character(ntaxa)
    current_block_type <- "num"
    in_block <- FALSE

    for (line in data_lines) {
      line <- trimws(line)
      if (line == "" || grepl("^;", line)) next

      if (grepl("^&", line)) {
        if (in_block) {
          if (block_taxon_idx != ntaxa + 1) warning("Block ended with incomplete taxa.")
          block_len <- nchar(current_block_seq[1])
          block_sizes <- c(block_sizes, block_len)
          block_types <- c(block_types, current_block_type)
          for (i in seq_len(ntaxa)) seq_accum[i] <- paste0(seq_accum[i], current_block_seq[i])
          current_block_seq <- character(ntaxa)
        }
        block_type_match <- regmatches(line, regexpr("(?i)(?<=&\\[)[^\\]]+", line, perl = TRUE))
        current_block_type <- if (length(block_type_match) > 0) block_type_match else "num"
        block_taxon_idx <- 1
        in_block <- TRUE
        remaining <- sub("^&\\[[^]]+\\]\\s*", "", line)
        if (nchar(remaining) > 0) {
          tokens <- unlist(strsplit(remaining, "\\s+"))
          for (i in seq(1, length(tokens), by = 2)) {
            if (i + 1 > length(tokens)) break
            taxon <- tokens[i]
            seq <- tokens[i + 1]
            if (block_taxon_idx <= ntaxa) {
              if (taxa_order[block_taxon_idx] == "") {
                taxa_order[block_taxon_idx] <- taxon
              } else if (taxa_order[block_taxon_idx] != taxon) warning("Taxon order mismatch")
              current_block_seq[block_taxon_idx] <- seq
              block_taxon_idx <- block_taxon_idx + 1
            }
          }
        }
        next
      }

      if (!in_block) {
        in_block <- TRUE
        current_block_type <- "num"
        block_taxon_idx <- 1
      }
      tokens <- unlist(strsplit(line, "\\s+"))
      for (i in seq(1, length(tokens), by = 2)) {
        if (i + 1 > length(tokens)) break
        taxon <- tokens[i]
        seq <- tokens[i + 1]
        if (block_taxon_idx <= ntaxa) {
          if (taxa_order[block_taxon_idx] == "") {
            taxa_order[block_taxon_idx] <- taxon
          } else if (taxa_order[block_taxon_idx] != taxon) warning("Taxon order mismatch")
          current_block_seq[block_taxon_idx] <- seq
          block_taxon_idx <- block_taxon_idx + 1
        }
      }
    }

    if (in_block && any(current_block_seq != "")) {
      block_len <- nchar(current_block_seq[1])
      block_sizes <- c(block_sizes, block_len)
      block_types <- c(block_types, current_block_type)
      for (i in seq_len(ntaxa)) seq_accum[i] <- paste0(seq_accum[i], current_block_seq[i])
    }

    total_len <- nchar(seq_accum[1])
    if (total_len != nchars) {
      warning("Total sequence length (", total_len, ") != header nchars (", nchars, "). Adjusting.")
      if (total_len < nchars) {
        pad <- paste(rep("?", nchars - total_len), collapse = "")
        for (i in seq_len(ntaxa)) seq_accum[i] <- paste0(seq_accum[i], pad)
        block_sizes <- c(block_sizes, nchars - total_len)
      } else {
        for (i in seq_len(ntaxa)) seq_accum[i] <- substr(seq_accum[i], 1, nchars)
        cum_len <- cumsum(block_sizes)
        keep <- cum_len <= nchars
        if (any(keep)) {
          block_sizes <- block_sizes[keep]
          if (sum(block_sizes) < nchars) block_sizes <- c(block_sizes, nchars - sum(block_sizes))
        } else {
          block_sizes <- nchars
        }
      }
    }

    token_pattern <- "\\[[^]]+\\]|\\{[^}]+\\}|\\([^)]*\\)|."
    mat_list <- regmatches(seq_accum, gregexpr(token_pattern, seq_accum, perl = TRUE))
    mat <- matrix(unlist(mat_list, use.names = FALSE), nrow = ntaxa, ncol = nchars, byrow = TRUE)
    rownames(mat) <- taxa_order

    attr(mat, "block_sizes") <- block_sizes
    attr(mat, "block_types") <- block_types
    mat
  }

  write_tnt_lines <- function(lines, filename) {
    if (length(lines) > 0 && lines[1] == "nstates 32") {
      text <- paste0(lines[1], "\n", paste(lines[-1], collapse = "\r\n"), "\r\n")
      writeBin(charToRaw(text), filename)
    } else {
      writeLines(lines, filename, sep = "\r\n")
    }
  }

  # ---------------------------------------------------------------------------
  # New helper: load a single input (file or R object) and return list with mat, type, ordered
  # ---------------------------------------------------------------------------
  load_single_input <- function(inp, fmt = NULL, name = "") {
    if (is.matrix(inp) || is.data.frame(inp)) {
      mat <- as.matrix(inp)
      type <- detect_matrix_type(mat)
      ordered <- integer(0)
    } else {
      if (is.null(fmt)) fmt <- detect_format_from_ext(inp)
      fmt <- match.arg(tolower(fmt), c("nexus", "tnt"))
      ordered <- get_ordered_chars(inp, fmt)
      if (fmt == "nexus") {
        header_type <- detect_type_from_nexus(inp)
        mat <- as.matrix(TreeTools::ReadCharacters(inp))
      } else {
        header_type <- detect_type_from_tnt(inp)
        mat <- read_tnt_manual(inp)
      }
      type <- if (header_type != "unknown") header_type else detect_matrix_type(mat)
    }
    list(mat = mat, type = type, ordered = ordered, name = name)
  }

  # ---------------------------------------------------------------------------
  # Determine list of inputs
  # ---------------------------------------------------------------------------
  if (!is.null(multiple_input)) {
    if (!dir.exists(multiple_input)) stop("Directory does not exist: ", multiple_input, call. = FALSE)
    files <- list.files(multiple_input, pattern = "\\.(tnt|nex|nexus)$", ignore.case = TRUE, full.names = TRUE)
    if (length(files) == 0) stop("No TNT or NEXUS files found in ", multiple_input)
    cat("Found", length(files), "files in", multiple_input, "\n")
    # Load each file
    inputs <- lapply(files, function(f) load_single_input(f, name = basename(f)))
  } else {
    # Traditional two-input mode
    if (is.null(input1) || is.null(input2)) stop("Both input1 and input2 must be provided (or use multiple_input).")
    inputs <- list(
      load_single_input(input1, input1_format, name = "matrix1"),
      load_single_input(input2, input2_format, name = "matrix2")
    )
  }

  # ---------------------------------------------------------------------------
  # Merge all matrices
  # ---------------------------------------------------------------------------
  all_mats <- lapply(inputs, `[[`, "mat")
  taxa_per_mat <- lapply(all_mats, rownames)
  all_taxa <- Reduce(union, taxa_per_mat)

  cat("\nMerging taxa across", length(inputs), "matrices\n")
  cat("Total unique taxa:", length(all_taxa), "\n")

  align_to_taxa <- function(mat, taxa) {
    out <- matrix("?",
      nrow = length(taxa), ncol = ncol(mat),
      dimnames = list(taxa, colnames(mat))
    )
    present <- intersect(taxa, rownames(mat))
    out[present, ] <- mat[present, , drop = FALSE]
    # Preserve block attributes
    attr(out, "block_sizes") <- attr(mat, "block_sizes")
    attr(out, "block_types") <- attr(mat, "block_types")
    out
  }

  aligned_mats <- lapply(all_mats, align_to_taxa, taxa = all_taxa)
  result <- do.call(cbind, aligned_mats)

  # Combine block attributes
  bs_all <- unlist(lapply(aligned_mats, function(m) {
    bs <- attr(m, "block_sizes")
    if (is.null(bs)) ncol(m) else bs
  }))
  bt_all <- unlist(lapply(seq_along(aligned_mats), function(i) {
    m <- aligned_mats[[i]]
    bt <- attr(m, "block_types")
    if (is.null(bt)) {
      rep(inputs[[i]]$type, length(attr(m, "block_sizes") %||% 1))
    } else {
      bt
    }
  }))
  attr(result, "block_sizes") <- bs_all
  attr(result, "block_types") <- bt_all

  # Combine ordered indices (NEXUS 1-based) with offsets
  combined_ordered <- integer(0)
  offset <- 0
  for (i in seq_along(inputs)) {
    ord <- inputs[[i]]$ordered
    if (length(ord) > 0) {
      combined_ordered <- c(combined_ordered, ord + offset)
    }
    offset <- offset + ncol(aligned_mats[[i]])
  }
  attr(result, "ordered") <- combined_ordered

  # For TNT writing we also need per-matrix ordered (0-based) with offsets
  tnt_ordered_combined <- integer(0)
  offset_tnt <- 0
  for (i in seq_along(inputs)) {
    ord <- inputs[[i]]$ordered
    if (length(ord) > 0) {
      tnt_ordered_combined <- c(tnt_ordered_combined, ord - 1 + offset_tnt)
    }
    offset_tnt <- offset_tnt + ncol(aligned_mats[[i]])
  }

  # Types for reporting
  types <- sapply(inputs, `[[`, "type")

  # ---------------------------------------------------------------------------
  # Print summary
  # ---------------------------------------------------------------------------
  cat("\nCharacters per matrix:\n")
  for (i in seq_along(aligned_mats)) {
    cat(sprintf("  %s: %d  [%s]", inputs[[i]]$name, ncol(aligned_mats[[i]]), toupper(types[i])))
    if (length(inputs[[i]]$ordered) > 0) cat(sprintf("  (%d ordered)", length(inputs[[i]]$ordered)))
    cat("\n")
  }
  cat("Total characters:", ncol(result))
  if (length(combined_ordered) > 0) cat(sprintf("  (including %d ordered/additive)", length(combined_ordered)))
  cat("\n\n")

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

      # Construct TNT output with all blocks
      ntaxa <- length(all_taxa)
      total_nchars <- ncol(result)
      lines <- c("xread", paste(total_nchars, ntaxa))

      write_tnt_block <- function(mat_part, bs, bt, start_col) {
        for (i in seq_along(bs)) {
          end <- start_col + bs[i] - 1
          submat <- mat_part[, start_col:end, drop = FALSE]
          header <- switch(bt[i],
            protein = "&[prot]",
            dna = "&[dna]",
            "&[num]"
          )
          lines <<- c(lines, "", header)
          for (j in seq_len(ntaxa)) {
            seq <- paste(format_tnt_sequence(submat[j, ]), collapse = "")
            lines <<- c(lines, paste0(rownames(mat_part)[j], "\t", seq))
          }
          start_col <- end + 1
        }
        start_col
      }

      col_cursor <- 1
      for (i in seq_along(aligned_mats)) {
        m <- aligned_mats[[i]]
        bs <- attr(m, "block_sizes")
        if (is.null(bs)) bs <- ncol(m)
        bt <- attr(m, "block_types")
        if (is.null(bt)) bt <- rep(types[i], length(bs))
        col_cursor <- write_tnt_block(result, bs, bt, col_cursor)
      }

      lines <- c(lines, ";", "")
      if (length(tnt_ordered_combined) > 0) {
        lines <- c(lines, paste0("cc + ", format_tnt_ordered(tnt_ordered_combined), " ;"), "")
      }
      lines <- c(lines, "proc /;", "comments 0", ";", "", "")
      write_tnt_lines(lines, output_file)

      cat("Concatenated matrix written to:", output_file, "\n")
      cat("Format: TNT")
    } else {
      if (!grepl("\\.(nexus|nex)$", output_file, ignore.case = TRUE)) output_file <- paste0(output_file, ".nexus")

      ape::write.nexus.data(result, file = output_file, format = "standard", interleaved = FALSE)
      tmp <- readLines(output_file)
      tmp <- gsub("INTERLEAVE=NO", "", tmp)
      tmp <- gsub("write.nexus.data.R", "RNODE", tmp)

      assump_block <- character(0)
      # Mixed if more than one unique type
      unique_types <- unique(types)
      if (length(unique_types) > 1) {
        # Build MIXED datatype string
        parts <- character(length(aligned_mats))
        start <- 1
        for (i in seq_along(aligned_mats)) {
          ncol_i <- ncol(aligned_mats[[i]])
          parts[i] <- paste0(toupper(types[i]), ":", start, "-", start + ncol_i - 1)
          start <- start + ncol_i
        }
        mixed_str <- paste0("DATATYPE=MIXED(", paste(parts, collapse = ", "), ")")
        tmp <- sub("DATATYPE=\\w+", mixed_str, tmp, ignore.case = TRUE)
      } else {
        tmp <- sub("DATATYPE=\\w+", paste0("DATATYPE=", toupper(unique_types)), tmp, ignore.case = TRUE)
      }

      if (length(combined_ordered) > 0) {
        assump_block <- c(
          assump_block,
          "BEGIN ASSUMPTIONS;",
          "    OPTIONS DEFTYPE=unord;",
          paste0("    TYPESET * Ordered = ", format_nexus_ordered(combined_ordered), ";"),
          "END;"
        )
      }
      if (length(assump_block) > 0) tmp <- c(tmp, assump_block)
      writeLines(tmp, output_file)

      cat("Concatenated matrix written to:", output_file, "\n")
      cat("Format: NEXUS")
    }
  }

  invisible(result)
}
