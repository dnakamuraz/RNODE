#' @title filterInvariants
#' @name filterInvariants
#' @description Delete invariant characters in a morphological matrix. We follow the definition of invariant site from IQ-Tree: (1) constant sites containing only a single character state in all sequences, (2) partially constant sites (N and/or -), and (3) ambiguously constant sites (e.g. C, Y and -).
#' @author Daniel YM Nakamura
#' @param input Input file (molecular or morphological matrix loaded locally or already loaded as a matrix object in R).
#' @param input_format To load from a local file: 'nexus' or 'tnt'.
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_onlyVARIANTS.nexus")
#' @examples
#' # Example
#' filterInvariants(input="../testdata/015_MORPH_data.nexus", output_index="../testdata/015_MORPH_data")
#'
#' @export
filterInvariants <- function(input, input_format, output_index) {

  # Load matrix
  if (is.matrix(input) || is.data.frame(input)) {
    # Already loaded matrix
    mat <- as.matrix(input)

  } else {
    # Read from a local file
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

  # --- Helper: parse states ---
  parse_state <- function(x) {
    x <- toupper(x)

    # Nucleotide ambiguity codes
    nuc_map <- list(
      A="A", C="C", G="G", T="T", U="T",
      R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
      K=c("G","T"), M=c("A","C"),
      B=c("C","G","T"), D=c("A","G","T"),
      H=c("A","C","T"), V=c("A","C","G"),
      N=c("A","C","G","T")
    )

    # Missing data
    if (x %in% c("?", "-")) return(NA)

    # Nucleotide?
    if (x %in% names(nuc_map)) return(nuc_map[[x]])

    # Morphological polymorphisms [01] {01} (01)
    if (grepl("^\\[[0-9]+\\]$", x) ||
        grepl("^\\{[0-9]+\\}$", x) ||
        grepl("^\\([0-9]+\\)$", x)) {
      inner <- gsub("\\[|\\]|\\{|\\}|\\(|\\)", "", x)
      return(strsplit(inner, "")[[1]])
    }

    # Single morphological digit
    if (grepl("^[0-9]$", x)) return(x)

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

    possible <- get_possible_states(col_no_na)

    # invariant only if intersection == 1 state
    !(length(possible) == 1)
  })

  deleted_columns <- which(!is_variable)

  # --- Report ---
  cat("Deleted", length(deleted_columns), "invariant columns.\n")
  if (length(deleted_columns) > 0)
    cat("Columns deleted:", paste(deleted_columns, collapse = ", "), "\n")

  # Filter to variable columns
  filtered <- mat[, is_variable, drop = FALSE]

  # --- Output file ---
  output_name <- paste0(output_index, "_onlyVARIANTS.nexus")

  # Keep the matrix orientation: rows are taxa, columns are characters.
  write.nexus.data(filtered, file = output_name, format = "standard", interleaved = FALSE)

  # Fix IQ-TREE incompatibility
  lines <- readLines(output_name)
  lines <- gsub("INTERLEAVE=NO", "", lines)
  writeLines(lines, output_name)

  return(filtered)
}
