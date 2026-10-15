#' @title consensusFAS
#' @name consensusFAS
#' @description Create a consensus sequence from DNA sequences in FASTA format.
#'   The consensus is computed column-wise using IUPAC ambiguity codes and appended
#'   to the alignment. Three consensus options are available:
#'   \itemize{
#'     \item \strong{relaxed}: The consensus nucleotide at each position is the
#'       IUPAC symbol that represents the union of all observed non-gap states.
#'       Gaps are ignored when other nucleotides are present; a column that is
#'       entirely gaps produces a gap in the consensus.
#'     \item \strong{conservative}: Same as relaxed, but positions where more than
#'       half of the sequences contain a gap are replaced with a question mark
#'       (\code{?}) in the consensus.
#'     \item \strong{both}: Appends both the relaxed and the conservative consensus
#'       sequences to the alignment.
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
#'   (default), \code{"conservative"}, or \code{"both"}. See Description
#'   for details.
#' @return A named character vector of aligned sequences with the consensus
#'   sequence(s) appended (named \code{Consensus}, \code{Consensus_relaxed},
#'   and/or \code{Consensus_conservative}).
#' @examples
#' # From a file, relaxed consensus
#' # result <- consensusFAS(input = "my_sequences.fasta")
#'
#' # Conservative consensus, write output
#' # result <- consensusFAS(input = "my_sequences.fasta",
#' #                        consensus = "conservative",
#' #                        output = "my_consensus.fasta")
#' @export
consensusFAS <- function(input, align = FALSE, output = NULL,
                         consensus = "relaxed") {

  # --- Validate consensus parameter ---
  consensus <- match.arg(consensus, choices = c("relaxed", "conservative", "both"))

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
  # 4. Build the consensus sequence(s)
  #    Both relaxed and conservative are computed in a single pass;
  #    only the requested one(s) are kept.
  # ---------------------------------------------------------------------------
  relaxed_chars      <- character(aln_length)
  conservative_chars <- character(aln_length)

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
      #   1. All-N                 -> N  (relaxed) / N  (conservative)
      #   2. Mix of N and gaps     -> N  (relaxed) / ?  (conservative)
      #   3. All gaps / ? / .      -> -  (relaxed) / -  (conservative)
      has_n   <- any(col == "N")
      has_gap <- any(col %in% c("-", ".", "?"))
      if (has_n && !has_gap) {
        # All-N
        relaxed_chars[j]      <- "N"
        conservative_chars[j] <- "N"
      } else if (has_n && has_gap) {
        # Mixed N and gaps
        relaxed_chars[j]      <- "N"
        conservative_chars[j] <- "?"
      } else {
        # All gaps / missing
        relaxed_chars[j]      <- "-"
        conservative_chars[j] <- "-"
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
  }

  relaxed_seq      <- paste(relaxed_chars, collapse = "")
  conservative_seq <- paste(conservative_chars, collapse = "")

  # ---------------------------------------------------------------------------
  # 5. Append consensus(es) to alignment
  # ---------------------------------------------------------------------------
  if (consensus == "both") {
    result <- c(seqs,
                Consensus_relaxed      = relaxed_seq,
                Consensus_conservative = conservative_seq)
    message("Consensus sequences (relaxed + conservative) computed successfully.")
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
