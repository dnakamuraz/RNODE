#' @title consensusFAS
#' @name consensusFAS
#' @description Create a consensus sequence from DNA sequences in FASTA format.
#'   The consensus is computed column-wise using IUPAC ambiguity codes and appended
#'   to the alignment. Four consensus options are available:
#'   \itemize{
#'     \item \strong{relaxed}: The consensus nucleotide at each position is the
#'       IUPAC symbol that represents the union of all observed non-gap states.
#'       Gaps are ignored when other nucleotides are present; a column that is
#'       entirely gaps produces a gap in the consensus.
#'     \item \strong{conservative}: Same as relaxed, but a position is replaced
#'       with \code{?} if more than half of the sequences are uninformative at
#'       that position (i.e., have \code{N}, \code{-}, or \code{?}). Note that
#'       columns where every sequence is \code{N} remain \code{N}, and columns
#'       with a mix of \code{N} and \code{-} (but no real nucleotide) become
#'       \code{?}.
#'     \item \strong{very_conservative}: Calls the nucleotide only when all
#'       informative sequences agree on a single base. If they disagree
#'       (ambiguity) but no gaps are present, inserts \code{N}; if they
#'       disagree and gaps are also present, inserts \code{?}.
#'     \item \strong{all}: Appends all three consensus sequences to the
#'       alignment (relaxed, conservative, and very_conservative).
#'   }
#' @author Daniel YM Nakamura
#' @param input A FASTA file path or a named character vector / \code{DNAbin}
#'   object already loaded in R. Each element should be one sequence.
#' @param align Logical; if \code{TRUE}, align the input sequences with MUSCLE
#'   (requires the \pkg{muscle} or system \code{muscle} executable). If
#'   \code{FALSE} (default), the sequences are assumed to be already aligned.
#' @param output If specified, the path to a local FASTA file to write the
#'   result (alignment + consensus). If \code{NULL} (default), the result is
#'   returned as an R object only.
#' @param consensus Method for building the consensus: \code{"relaxed"}
#'   (default), \code{"conservative"}, \code{"very_conservative"}, or
#'   \code{"all"}. See Description for details.
#' @param trimming Logical; if \code{TRUE}, flanking columns where every
#'   sequence contains only \code{N}, \code{-}, or \code{?} are removed from
#'   both the alignment and the consensus before output. Internal uninformative
#'   columns are not affected. Default is \code{FALSE}.
#' @return A named character vector of aligned sequences with the consensus
#'   sequence(s) appended (named \code{Consensus}, \code{Consensus_relaxed},
#'   \code{Consensus_conservative}, and/or \code{Consensus_very_conservative}).
#' @examples
#' # From a file, relaxed consensus
#' # result <- consensusFAS(input = "my_sequences.fasta")
#'
#' # Conservative consensus, written to file, with end-trimming
#' # result <- consensusFAS(input = "my_sequences.fasta",
#' #                        consensus = "conservative",
#' #                        trimming = TRUE,
#' #                        output = "my_consensus.fasta")
#' @export
consensusFAS <- function(input, align = FALSE, output = NULL,
                         consensus = "relaxed", trimming = FALSE) {

  # --- Validate consensus parameter ---
  consensus <- match.arg(consensus, choices = c("relaxed", "conservative",
                                                 "very_conservative", "all"))

  # ---------------------------------------------------------------------------
  # IUPAC lookup: sorted nucleotide set -> IUPAC code
  # ---------------------------------------------------------------------------
  iupac_map <- list(
    "A"    = "A",
    "C"    = "C",
    "G"    = "G",
    "T"    = "T",
    "AG"   = "R",
    "CT"   = "Y",
    "CG"   = "S",
    "AT"   = "W",
    "GT"   = "K",
    "AC"   = "M",
    "CGT"  = "B",
    "AGT"  = "D",
    "ACT"  = "H",
    "ACG"  = "V",
    "ACGT" = "N"
  )

  # Reverse IUPAC map: symbol -> constituent nucleotides
  iupac_to_bases <- list(
    A = "A", C = "C", G = "G", T = "T", U = "T",
    R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"), W = c("A", "T"),
    K = c("G", "T"), M = c("A", "C"),
    B = c("C", "G", "T"), D = c("A", "G", "T"),
    H = c("A", "C", "T"), V = c("A", "C", "G"),
    N = c("A", "C", "G", "T")
  )

  # ---------------------------------------------------------------------------
  # 1. Read input sequences
  # ---------------------------------------------------------------------------
  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    # Read FASTA from file
    seqs <- read_fasta(input)
  } else if (inherits(input, "DNAbin")) {
    # Convert ape::DNAbin to named character vector
    if (is.matrix(input)) {
      seqs <- apply(input, 1, function(row) {
        paste(as.character(row), collapse = "")
      })
    } else if (is.list(input)) {
      seqs <- sapply(input, function(s) paste(as.character(s), collapse = ""))
    } else {
      stop("Unsupported DNAbin format.")
    }
    seqs <- toupper(seqs)
  } else if (is.character(input) && !is.null(names(input))) {
    # Already a named character vector
    seqs <- toupper(input)
  } else {
    stop("'input' must be a path to a FASTA file, a named character vector, or a DNAbin object.")
  }

  if (length(seqs) < 2) {
    stop("At least two sequences are required to compute a consensus.")
  }

  # ---------------------------------------------------------------------------
  # 2. Optionally align sequences
  # ---------------------------------------------------------------------------
  if (align) {
    message("Aligning sequences with MUSCLE...")

    # Write sequences to a temporary FASTA file
    tmp_in  <- tempfile(fileext = ".fasta")
    tmp_out <- tempfile(fileext = ".fasta")
    write_fasta(seqs, tmp_in)

    # Try system MUSCLE first
    muscle_cmd <- Sys.which("muscle")
    if (nzchar(muscle_cmd)) {
      # Detect MUSCLE version (v5+ uses -align/-output, v3 uses -in/-out)
      version_check <- tryCatch(
        system2(muscle_cmd, args = "-version", stdout = TRUE, stderr = TRUE),
        error = function(e) ""
      )
      version_str <- paste(version_check, collapse = " ")

      if (grepl("muscle\\s+[5-9]", version_str, ignore.case = TRUE) ||
          grepl("muscle\\s+v[5-9]", version_str, ignore.case = TRUE)) {
        # MUSCLE v5+
        ret <- system2(muscle_cmd,
                       args = c("-align", tmp_in, "-output", tmp_out),
                       stdout = TRUE, stderr = TRUE)
      } else {
        # MUSCLE v3 / other
        ret <- system2(muscle_cmd,
                       args = c("-in", tmp_in, "-out", tmp_out),
                       stdout = TRUE, stderr = TRUE)
      }

      if (!file.exists(tmp_out) || file.size(tmp_out) == 0) {
        stop("MUSCLE alignment failed. Ensure MUSCLE is installed and in your PATH.\n",
             paste(ret, collapse = "\n"))
      }

      seqs <- read_fasta(tmp_out)
    } else if (requireNamespace("muscle", quietly = TRUE)) {
      # Use the R muscle package as fallback
      dna <- Biostrings::DNAStringSet(seqs)
      aln <- muscle::muscle(dna)
      seqs <- setNames(as.character(aln), names(seqs))
    } else {
      stop("MUSCLE is not available. Install the system 'muscle' executable ",
           "or the R 'muscle' package.")
    }

    unlink(c(tmp_in, tmp_out))
    message("Alignment complete.")
  }

  # ---------------------------------------------------------------------------
  # 3. Validate alignment (equal length)
  # ---------------------------------------------------------------------------
  seq_lengths <- nchar(seqs)
  if (length(unique(seq_lengths)) != 1) {
    stop("Sequences are not all the same length. ",
         "Set align = TRUE to align them, or provide a pre-aligned FASTA file.")
  }

  aln_length <- seq_lengths[1]
  n_seqs     <- length(seqs)

  # Split sequences into a matrix of characters
  seq_mat <- do.call(rbind, strsplit(seqs, ""))

  # ---------------------------------------------------------------------------
  # 4. Optional: trim flanking uninformative columns
  #    A column is "uninformative" if every sequence has N, -, or ? there.
  #    Only leading and trailing runs of such columns are removed.
  # ---------------------------------------------------------------------------
  if (trimming) {
    uninformative_col <- function(mat_col) {
      all(toupper(mat_col) %in% c("N", "-", "?", "."))
    }

    col_is_uninf <- sapply(seq_len(aln_length), function(j) {
      uninformative_col(seq_mat[, j])
    })

    # Find first and last informative column
    first_inf <- which(!col_is_uninf)[1]
    last_inf  <- tail(which(!col_is_uninf), 1)

    if (is.na(first_inf)) {
      warning("All columns are uninformative after trimming; returning untrimmed alignment.")
    } else {
      keep_cols  <- first_inf:last_inf
      seq_mat    <- seq_mat[, keep_cols, drop = FALSE]
      aln_length <- length(keep_cols)
      # Rebuild sequence strings to match trimmed matrix
      seqs <- apply(seq_mat, 1, paste, collapse = "")
      message("Trimming: kept columns ", first_inf, "-", last_inf,
              " (", aln_length, " bp after trimming).")
    }
  }

  # ---------------------------------------------------------------------------
  # 5. Build the consensus sequence(s)
  #    Both relaxed and conservative are computed in a single pass;
  #    only the requested one(s) are kept.
  # ---------------------------------------------------------------------------
  relaxed_chars          <- character(aln_length)
  conservative_chars     <- character(aln_length)
  very_conservative_chars <- character(aln_length)

  for (j in seq_len(aln_length)) {
    col <- toupper(seq_mat[, j])

    # Expand each character to its constituent bases (handle IUPAC in input)
    # N is treated as missing data (unknown base), not as {A,C,G,T}
    base_sets <- lapply(col, function(ch) {
      if (ch %in% c("-", ".", "?", "N")) {
        return(NULL)  # treat as gap / missing / unknown
      }
      up <- toupper(ch)
      if (!is.null(iupac_to_bases[[up]])) {
        return(iupac_to_bases[[up]])
      }
      return(NULL)  # unknown characters treated as missing
    })

    # Count how many sequences have a non-gap character
    non_gap_idx <- which(!sapply(base_sets, is.null))
    n_non_gap   <- length(non_gap_idx)

    if (n_non_gap == 0) {
      # Three sub-cases when no informative base is present:
      #   1. All-N                 -> N  (relaxed) / N  (conservative) / N  (very_conservative)
      #   2. Mix of N and gaps     -> N  (relaxed) / ?  (conservative) / ?  (very_conservative)
      #   3. All gaps / ? / .      -> -  (relaxed) / -  (conservative) / -  (very_conservative)
      has_n   <- any(col == "N")
      has_gap <- any(col %in% c("-", ".", "?"))
      if (has_n && !has_gap) {
        # All-N
        relaxed_chars[j]           <- "N"
        conservative_chars[j]      <- "N"
        very_conservative_chars[j] <- "N"
      } else if (has_n && has_gap) {
        # Mixed N and gaps
        relaxed_chars[j]           <- "N"
        conservative_chars[j]      <- "?"
        very_conservative_chars[j] <- "?"
      } else {
        # All gaps / missing
        relaxed_chars[j]           <- "-"
        conservative_chars[j]      <- "-"
        very_conservative_chars[j] <- "-"
      }
      next
    }

    # Union of all observed bases
    all_bases <- sort(unique(unlist(base_sets[non_gap_idx])))
    key       <- paste(all_bases, collapse = "")

    # Look up the IUPAC symbol
    iupac_symbol <- iupac_map[[key]]
    if (is.null(iupac_symbol)) {
      iupac_symbol <- "N"  # fallback for unexpected combinations
    }

    # Relaxed: always use the IUPAC symbol
    relaxed_chars[j] <- iupac_symbol

    # Conservative: use "?" if more than half the sequences are gaps/missing
    n_gaps <- n_seqs - n_non_gap
    if (n_gaps > n_seqs / 2) {
      conservative_chars[j] <- "?"
    } else {
      conservative_chars[j] <- iupac_symbol
    }

    # Very conservative: only call the base when all informative seqs agree
    #   (i.e., the IUPAC symbol is a single unambiguous base).
    #   If ambiguous AND gaps present -> "?"; if ambiguous but no gaps -> "N".
    is_ambiguous <- length(all_bases) > 1
    has_gap_col  <- n_non_gap < n_seqs
    if (!is_ambiguous) {
      very_conservative_chars[j] <- iupac_symbol  # single base, unambiguous
    } else if (has_gap_col) {
      very_conservative_chars[j] <- "?"  # disagreement + gaps
    } else {
      very_conservative_chars[j] <- "N"  # disagreement, no gaps
    }
  }

  relaxed_seq           <- paste(relaxed_chars, collapse = "")
  conservative_seq      <- paste(conservative_chars, collapse = "")
  very_conservative_seq <- paste(very_conservative_chars, collapse = "")

  # ---------------------------------------------------------------------------
  # 6. Append consensus(es) to alignment
  # ---------------------------------------------------------------------------
  if (consensus == "all") {
    result <- c(seqs,
                Consensus_relaxed           = relaxed_seq,
                Consensus_conservative      = conservative_seq,
                Consensus_very_conservative = very_conservative_seq)
    message("Consensus sequences (relaxed + conservative + very_conservative) computed successfully.")
  } else if (consensus == "very_conservative") {
    result <- c(seqs, Consensus = very_conservative_seq)
    message("Consensus sequence (very_conservative) computed successfully.")
  } else if (consensus == "conservative") {
    result <- c(seqs, Consensus = conservative_seq)
    message("Consensus sequence (conservative) computed successfully.")
  } else {
    result <- c(seqs, Consensus = relaxed_seq)
    message("Consensus sequence (relaxed) computed successfully.")
  }

  message("Consensus length: ", aln_length, " bp")

  # ---------------------------------------------------------------------------
  # 6. Optionally write output
  # ---------------------------------------------------------------------------
  if (!is.null(output)) {
    write_fasta(result, output)
    message("Output written to: ", output)
  }

  return(invisible(result))
}


# =============================================================================
# Internal helpers: read / write FASTA
# =============================================================================

#' Read a FASTA file into a named character vector
#' @param filepath Path to the FASTA file.
#' @return A named character vector (names = sequence headers, values = sequences).
#' @keywords internal
read_fasta <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  # Remove empty lines
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0) {
    stop("The FASTA file is empty: ", filepath)
  }

  names_vec <- character()
  seqs_vec  <- character()
  current_name <- NULL
  current_seq  <- ""

  for (line in lines) {
    if (startsWith(line, ">")) {
      # Save previous sequence
      if (!is.null(current_name)) {
        names_vec  <- c(names_vec, current_name)
        seqs_vec   <- c(seqs_vec, current_seq)
      }
      current_name <- sub("^>\\s*", "", line)
      current_seq  <- ""
    } else {
      current_seq <- paste0(current_seq, trimws(line))
    }
  }
  # Save last sequence
  if (!is.null(current_name)) {
    names_vec <- c(names_vec, current_name)
    seqs_vec  <- c(seqs_vec, current_seq)
  }

  seqs_vec <- toupper(seqs_vec)
  names(seqs_vec) <- names_vec
  return(seqs_vec)
}


#' Write a named character vector to a FASTA file
#' @param seqs Named character vector of sequences.
#' @param filepath Path to the output FASTA file.
#' @keywords internal
write_fasta <- function(seqs, filepath) {
  out_lines <- character()
  for (i in seq_along(seqs)) {
    out_lines <- c(out_lines,
                   paste0(">", names(seqs)[i]),
                   seqs[i])
  }
  writeLines(out_lines, filepath)
}
