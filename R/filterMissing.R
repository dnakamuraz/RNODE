#' @title filterMissing
#' @name filterMissing
#' @description Filter taxa (rows) and/or characters (columns) in a matrix containing only missing data.
#' @author Daniel YM Nakamura
#' @param input Input file (morphological matrix in Nexus format).
#' @param input_format To load from a local file: 'nexus' or 'tnt'.
#' @param output_path Output path (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @param output_format Format to write the output matrix. Options: "nexus" (default) or "tnt".
#' @param missing Parameter specifying if rows and/or columns in which all cells are missing data (?) should be removed. Options: "row" (default i.e. terminals), "column" (i.e. characters, transformation series), "both".
#' @examples
#' # Example
#' filterMissing (input="testdata/test_filterMissing.nexus", output_path="testdata/test", missing="row")
#'
#' @export
filterMissing = function(input, input_format = NULL,
                       output_path,
                       output_format="nexus",
                       missing="row") {

  # Load matrix
  if (is.matrix(input) || is.data.frame(input)) {
    # Already loaded matrix
    data <- as.matrix(input)

  } else {
    # Read from a local file
    if (is.null(input_format)) {
      stop("If 'input' is a file path, you must provide input_format = 'nexus' or 'tnt'")
    }
    if (input_format == "nexus") {
      data <- TreeTools::ReadCharacters(input)
    } else if (input_format == "tnt") {
      # Helper function to read TNT file manually preserving exact format
      read_tnt_manual <- function(tnt_file) {
        lines <- readLines(tnt_file)
        
        # Remove preamble lines (like "nstates prot;")
        xread_idx <- which(grepl("^xread", lines, ignore.case = TRUE))[1]
        if (is.na(xread_idx)) {
          stop("Could not find 'xread' block in TNT file: ", tnt_file)
        }
        lines <- lines[xread_idx:length(lines)]
        
        # Parse header: line 1 is xread, line 2 is nchars ntaxa
        header_parts <- strsplit(trimws(lines[2]), "\\s+")[[1]]
        nchars <- as.numeric(header_parts[1])
        ntaxa <- as.numeric(header_parts[2])
        if (is.na(nchars) || is.na(ntaxa)) {
          stop("Could not parse TNT dimensions in file: ", tnt_file)
        }
        
        # Keep only taxon sequence lines. TNT block markers such as "& [num]"
        # and the trailing command section are not part of the matrix.
        data_lines <- trimws(lines[-c(1, 2)])
        data_lines <- data_lines[nzchar(data_lines)]
        data_lines <- data_lines[!grepl("^&\\s*\\[", data_lines)]
        end_idx <- which(data_lines == ";")[1]
        if (!is.na(end_idx)) {
          data_lines <- data_lines[seq_len(end_idx - 1)]
        }
        
        if (length(data_lines) < ntaxa) {
          stop("Expected ", ntaxa, " taxa, but found ", length(data_lines),
               " sequence lines in TNT file: ", tnt_file)
        }
        seq_lines <- data_lines[seq_len(ntaxa)]
        
        seq_parts <- regmatches(
          seq_lines,
          regexec("^([^\\s]+)\\s+(.+)$", seq_lines, perl = TRUE)
        )
        parsed <- vapply(seq_parts, length, integer(1)) == 3
        if (!all(parsed)) {
          stop("Could not parse taxon and sequence on line(s): ",
               paste(which(!parsed), collapse = ", "))
        }
        
        taxa_names <- vapply(seq_parts, `[`, character(1), 2)
        sequences <- vapply(seq_parts, `[`, character(1), 3)
        sequences <- gsub("\\s+", "", sequences, perl = TRUE)
        
        # Fast path for simple aligned sequences
        if (!any(grepl("[\\[\\{\\(]", sequences, perl = TRUE))) {
          seq_widths <- nchar(sequences, type = "chars")
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
          if (parsed_nchars != nchars) {
            warning("TNT header declares ", nchars, " characters, but parsed ",
                    parsed_nchars, " characters in file: ", tnt_file)
          }
          mat <- matrix(
            unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE),
            nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE
          )
          rownames(mat) <- taxa_names
          return(mat)
        }
        
        # Slower path for polymorphic/inapplicable tokens
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
        mat <- matrix(unlist(mat_list, use.names = FALSE),
                      nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
        rownames(mat) <- taxa_names
        return(mat)
      }
      
      data <- read_tnt_manual(input)
    } else {
      stop("input_format must be 'nexus' or 'tnt'")
    }
    data <- as.matrix(data)
  }

  # If missing is 'row', delete rows containing only ?
  if (missing == 'row' || missing == "both") {
    rows_to_delete <- apply(data, 1, function(row) all(row == "?"))
    cat("Rows deleted:", which(rows_to_delete), "\n")
    data <- data[!rows_to_delete, ]
  }

  # If missing is 'row', delete rows containing only ?
  if (missing == 'column' || missing == "both") {
    cols_to_delete <- apply(data, 2, function(col) all(col == "?"))
    cat("Columns deleted:", which(cols_to_delete), "\n")
    data <- data[!cols_to_delete, ]
  }

  # Validate output_format
  if (!is.null(output_format)) {
    output_format <- tolower(output_format)
  }
  
  if (!is.null(output_format) && output_format == "tnt") {
    name = paste0(output_path, "_FILTERED.tnt")
    # TNT uses [01] for ambiguities; convert any (01)-style tokens before writing
    data[] <- gsub("\\(([^)]+)\\)", "[\\1]", data, perl = TRUE)
    TreeTools::WriteTntCharacters(data, filepath = name)
  } else {
    # Default to nexus
    name = paste0(output_path, "_FILTERED.nexus")
    # Write a temporary file
    ape::write.nexus.data(data, file=name, format="standard", interleaved=F)
    # Remove the interleave section to avoid errors in IQTREE
    temp = gsub("INTERLEAVE=NO", "", readLines(name))
    # Write a permanent file
    writeLines(temp, name)
  }
}
