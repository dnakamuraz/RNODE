#' @title filterInvariants
#' @name filterInvariants
#' @description Delete invariant characters in a morphological matrix. We follow the definition of invariant site from IQ-Tree: (1) constant sites containing only a single character state in all sequences, (2) partially constant sites (N and/or -), and (3) ambiguously constant sites (e.g. C, Y and -).
#' @author Daniel YM Nakamura
#' @param input Input file (molecular or morphological matrix loaded locally or already loaded as a matrix object in R).
#' @param input_format To load from a local file: 'nexus', 'tnt', or NULL for auto-detect. Default is NULL.
#' @param output_format To write a local file: 'tnt', 'nexus', or NULL (default; only returns as an R object).
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_onlyVARIANTS.nexus"). Required if output_format is not NULL.
#' @examples
#' # Example returning only R object:
#' mat_filtered <- filterInvariants(input = "../testdata/015_MORPH_data.nexus")
#'
#' # Example writing locally as TNT:
#' filterInvariants(input = "../testdata/015_MORPH_data.nexus", output_format = "tnt", output_index = "../testdata/015_MORPH_data")
#'
#' @export
filterInvariants <- function(input, input_format = NULL, output_format = NULL, output_index = NULL) {
  # Load matrix
  if (is.matrix(input) || is.data.frame(input)) {
    # Already loaded matrix
    mat <- as.matrix(input)
  } else {
    # Read from a local file

    # --- Auto-detect format if input_format is NULL ---
    if (is.null(input_format)) {
      # Read the first 20 lines to guess the format
      first_lines <- readLines(input, n = 20, warn = FALSE)
      # Clean whitespaces
      first_lines <- trimws(first_lines)
      first_lines <- first_lines[first_lines != ""]

      if (length(first_lines) == 0) {
        stop("The input file appears to be empty.")
      }

      # 1. Check for NEXUS
      if (any(grepl("(?i)^#NEXUS", first_lines))) {
        input_format <- "nexus"
        message("Auto-detected input format: NEXUS")

        # 2. Check for TNT keywords
      } else if (any(grepl("(?i)^(xread|nstates|mxram|rseed)", first_lines))) {
        input_format <- "tnt"
        message("Auto-detected input format: TNT")

        # 3. Fallback to file extension
      } else {
        if (grepl("\\.nex(us)?$", input, ignore.case = TRUE)) {
          input_format <- "nexus"
          message("Auto-detected input format from extension: NEXUS")
        } else if (grepl("\\.tnt$", input, ignore.case = TRUE)) {
          input_format <- "tnt"
          message("Auto-detected input format from extension: TNT")
        } else {
          stop("Could not automatically determine input_format. Please explicitly provide input_format = 'nexus' or 'tnt'")
        }
      }
    }
    # --------------------------------------------------

    if (input_format == "nexus") {
      mat <- TreeTools::ReadCharacters(input)
    } else if (input_format == "tnt") {
      # Read lines to strip TNT block tags like &[num], &[dna], &[prot], etc.
      # TreeTools natively can misinterpret these tags as character states shifting
      # the alignment, leading to 'Unrecognized symbol: [NUM]' errors.
      lines <- readLines(input, warn = FALSE)
      cleaned_lines <- gsub("&\\s*\\[\\s*[a-zA-Z]+\\s*\\]", "", lines)

      temp_input <- tempfile(fileext = ".tnt")
      writeLines(cleaned_lines, temp_input)

      mat <- TreeTools::ReadTntCharacters(temp_input)
      unlink(temp_input) # remove temp file
    } else {
      stop("input_format must be 'nexus' or 'tnt'")
    }
    mat <- as.matrix(mat)
  }

  # --- Helper: parse states ---
  parse_state <- function(x) {
    x <- toupper(x)

    # Nucleotide ambiguity codes
    nuc_map <- list(
      A = "A", C = "C", G = "G", T = "T", U = "T",
      R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
      K = c("G", "T"), M = c("A", "C"),
      B = c("C", "G", "T"), D = c("A", "G", "T"),
      H = c("A", "C", "T"), V = c("A", "C", "G"),
      N = c("A", "C", "G", "T")
    )

    # Missing data
    if (x %in% c("?", "-")) {
      return(NA)
    }

    # Nucleotide?
    if (x %in% names(nuc_map)) {
      return(nuc_map[[x]])
    }

    # Polymorphisms [01] {01} (01) - Extended to A-Z for >9 state morph matrices & Nucleotide polymorphisms
    if (grepl("^\\[[0-9A-Z]+\\]$", x) ||
      grepl("^\\{[0-9A-Z]+\\}$", x) ||
      grepl("^\\([0-9A-Z]+\\)$", x)) {
      inner <- gsub("\\[|\\]|\\{|\\}|\\(|\\)", "", x)
      return(strsplit(inner, "")[[1]])
    }

    # Single morphological digit or character (Extended for >9 states where Letters are used)
    if (grepl("^[0-9A-Z]$", x)) {
      return(x)
    }

    stop(paste("Unrecognized symbol:", x))
  }

  morph_states <- as.character(0:9)

  # Determine possible states in a column
  get_possible_states <- function(chars) {
    sets <- lapply(chars, parse_state)

    # Missing data (NA) expands to all morph states
    sets <- lapply(sets, function(s) {
      if (length(s) == 1 && is.na(s)) morph_states else s
    })

    Reduce(intersect, sets)
  }

  # Determine variable vs invariant columns
  is_variable <- sapply(seq_len(ncol(mat)), function(i) {
    col <- mat[, i]
    col_no_na <- col[!is.na(col)]

    # Prevent failure if the entire column contains only NAs
    if (length(col_no_na) == 0) {
      return(FALSE) # Fully missing columns are invariant
    }

    possible <- get_possible_states(col_no_na)

    # invariant only if intersection == 1 state
    !(length(possible) == 1)
  })

  deleted_columns <- which(!is_variable)

  # --- Report ---
  cat("Deleted", length(deleted_columns), "invariant columns.\n")
  if (length(deleted_columns) > 0) {
    cat("Columns deleted:", paste(deleted_columns, collapse = ", "), "\n")
  }

  # Filter to variable columns
  filtered <- mat[, is_variable, drop = FALSE]

  # --- Output file (If requested) ---
  if (!is.null(output_format)) {
    if (is.null(output_index)) {
      stop("You must provide an 'output_index' when 'output_format' is specified.")
    }

    if (ncol(filtered) == 0) {
      warning("All columns were evaluated as invariant. Skipping writing local file to avoid errors.")
      return(filtered)
    }

    fmt <- tolower(output_format)

    if (fmt == "nexus") {
      output_name <- paste0(output_index, "_onlyVARIANTS.nexus")
      ape::write.nexus.data(filtered, file = output_name, format = "standard", interleaved = FALSE)

      # Fix IQ-TREE incompatibility
      lines_out <- readLines(output_name)
      lines_out <- gsub("INTERLEAVE=NO", "", lines_out)
      writeLines(lines_out, output_name)

      message("Saved NEXUS output to: ", output_name)
    } else if (fmt == "tnt") {
      output_name <- paste0(output_index, "_onlyVARIANTS.tnt")

      # Prepare data for TNT
      n_taxa <- nrow(filtered)
      n_chars <- ncol(filtered)

      # Translate NA back to standard TNT ? code
      filtered_out <- filtered
      filtered_out[is.na(filtered_out)] <- "?"

      taxa_names <- rownames(filtered_out)
      if (is.null(taxa_names)) taxa_names <- paste0("Taxon", seq_len(n_taxa))

      # Calculate padding for clean alignment
      pad <- max(nchar(taxa_names)) + 2

      # Write out standard TNT
      cat("xread\n", file = output_name)
      cat("'Generated by filterInvariants'\n", file = output_name, append = TRUE)
      cat(sprintf("%d %d\n", n_chars, n_taxa), file = output_name, append = TRUE)

      for (i in seq_len(n_taxa)) {
        seq_string <- paste0(filtered_out[i, ], collapse = "")
        cat(sprintf("%-*s %s\n", pad, taxa_names[i], seq_string), file = output_name, append = TRUE)
      }
      cat(";\n", file = output_name, append = TRUE)

      message("Saved TNT output to: ", output_name)
    } else {
      stop("output_format must be 'nexus', 'tnt', or NULL.")
    }
  }

  return(filtered)
}
