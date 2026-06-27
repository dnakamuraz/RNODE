#' @title filterMissing
#' @name filterMissing
#' @description Filter taxa (rows) and/or characters (columns) in a matrix containing only missing data.
#' @author Daniel YM Nakamura
#' @param input Input file (morphological matrix in Nexus format).
#' @param input_format To load from a local file: 'nexus' or 'tnt'.
#' @param output Output path (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @param output_format Format to write the output matrix. Options: "nexus" (default) or "tnt".
#' @param missing Parameter specifying if rows and/or columns in which all cells are missing data (?) should be removed. Options: "row" (default i.e. terminals), "column" (i.e. characters, transformation series), "both".
#' @examples
#' \dontrun{
#' # Example
#' filterMissing(input = "testdata/test_filterMissing.nexus", output = "testdata/test", missing = "row")
#' }
#'
#' @export
filterMissing <- function(input, input_format = NULL, output, output_format = "nexus", missing = "row") {
  # Load matrix
  if (is.matrix(input) || is.data.frame(input)) {
    # Already loaded matrix
    data <- as.matrix(input)
  } else if (is.character(input)) {
    # Ensure input path actually exists
    if (!file.exists(input)) {
      stop("The input file path does not exist.")
    }

    # Auto-detect format if set to NULL
    if (is.null(input_format)) {
      ext <- tolower(tools::file_ext(input))
      if (ext %in% c("nex", "nexus")) {
        input_format <- "nexus"
      } else if (ext %in% c("tnt")) {
        input_format <- "tnt"
      } else {
        peek_lines <- readLines(input, n = 20, warn = FALSE)
        if (any(grepl("^#NEXUS", peek_lines, ignore.case = TRUE))) {
          input_format <- "nexus"
        } else if (any(grepl("^(xread|nstates)", peek_lines, ignore.case = TRUE))) {
          input_format <- "tnt"
        } else {
          stop("Could not automatically detect the input format from the file. Please specify input_format = 'nexus' or 'tnt'.")
        }
      }
    }

    if (input_format == "nexus") {
      data <- TreeTools::ReadCharacters(input)
    } else if (input_format == "tnt") {
      # Helper function to read TNT file manually and preserve original structure
      read_tnt_manual <- function(tnt_file) {
        raw_lines <- readLines(tnt_file)

        xread_idx <- which(grepl("^xread", raw_lines, ignore.case = TRUE))[1]
        if (is.na(xread_idx)) stop("Could not find 'xread' block in TNT file: ", tnt_file)

        # Save preamble (e.g. nstates, comments before matrix)
        preamble <- if (xread_idx > 1) raw_lines[seq_len(xread_idx - 1)] else character(0)

        lines <- raw_lines[xread_idx:length(raw_lines)]

        header_parts <- strsplit(trimws(lines[2]), "\\s+")[[1]]
        nchars <- as.numeric(header_parts[1])
        ntaxa <- as.numeric(header_parts[2])
        if (is.na(nchars) || is.na(ntaxa)) {
          stop("Could not parse TNT dimensions in file: ", tnt_file)
        }

        # Find closing semicolon of the matrix block
        semi_matches <- which(trimws(lines) == ";")
        matrix_semi_idx <- semi_matches[semi_matches > 2][1]

        # Extract everything after the main matrix as the trailer (cc blocks, proc /;)
        if (!is.na(matrix_semi_idx) && matrix_semi_idx < length(lines)) {
          trailer <- lines[(matrix_semi_idx + 1):length(lines)]
        } else {
          trailer <- character(0)
        }

        # Look for tag blocks like &[num] before the data
        matrix_body <- trimws(lines[3:(matrix_semi_idx - 1)])
        tag_match <- grep("^&\\s*\\[", matrix_body, value = TRUE)
        block_tag <- if (length(tag_match) > 0) tag_match[1] else character(0)

        # Filter for actual sequences
        seq_lines <- matrix_body[nzchar(matrix_body)]
        seq_lines <- seq_lines[!grepl("^&\\s*\\[", seq_lines)]

        if (length(seq_lines) < ntaxa) stop("Expected ", ntaxa, " sequence lines.")
        seq_lines <- seq_lines[seq_len(ntaxa)]

        seq_parts <- regmatches(seq_lines, regexec("^([^\\s]+)\\s+(.+)$", seq_lines, perl = TRUE))
        parsed <- vapply(seq_parts, length, integer(1)) == 3
        if (!all(parsed)) stop("Could not parse taxon and sequence on line(s): ", paste(which(!parsed), collapse = ", "))

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
              nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE
            )
            rownames(mat) <- taxa_names
          } else {
            mat <- matrix(unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE),
              nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE
            )
            rownames(mat) <- taxa_names
          }
        } else {
          # Slower path for polymorphic tokens
          token_pattern <- "\\[[^]]+\\]|\\{[^}]+\\}|\\([^)]*\\)|."
          mat_list <- regmatches(sequences, gregexpr(token_pattern, sequences, perl = TRUE))
          token_counts <- lengths(mat_list)
          parsed_nchars <- max(token_counts)
          if (!all(token_counts == parsed_nchars)) {
            mat_list <- lapply(mat_list, function(x) {
              if (length(x) < parsed_nchars) c(x, rep("?", parsed_nchars - length(x))) else x
            })
          }
          mat <- matrix(unlist(mat_list, use.names = FALSE), nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
          rownames(mat) <- taxa_names
        }

        # Attach the extracted format pieces as attributes
        attr(mat, "tnt_preamble") <- preamble
        attr(mat, "tnt_trailer") <- trailer
        attr(mat, "tnt_block_tag") <- block_tag
        return(mat)
      }

      data <- read_tnt_manual(input)
    } else {
      stop("input_format must be 'nexus' or 'tnt'")
    }
    data <- as.matrix(data)
  } else {
    stop("Input must be an R matrix/data frame or a valid file path character string.")
  }

  # Save attributes before subsetting (R strips attributes when matrices are subset)
  tnt_preamble <- attr(data, "tnt_preamble")
  tnt_trailer <- attr(data, "tnt_trailer")
  tnt_block_tag <- attr(data, "tnt_block_tag")

  cols_to_delete <- logical(ncol(data))

  # Delete rows containing only '?'
  if (missing == "row" || missing == "both") {
    rows_to_delete <- apply(data, 1, function(row) all(row == "?"))
    if (any(rows_to_delete)) {
      cat("Rows deleted:", which(rows_to_delete), "\n")
      data <- data[!rows_to_delete, , drop = FALSE]
    }
  }

  # Delete columns containing only '?'
  if (missing == "column" || missing == "both") {
    cols_to_delete <- apply(data, 2, function(col) all(col == "?"))
    if (any(cols_to_delete)) {
      cat("Columns deleted:", which(cols_to_delete), "\n")
      data <- data[, !cols_to_delete, drop = FALSE]
    }
  }

  # Warn if deleting columns breaks TNT trailing blocks or ordered character assumptions
  # NOTE: Added dual backward slashes (\\) so grepl's regex successfully catches spaced characters
  if (any(cols_to_delete) && !is.null(tnt_trailer) && any(grepl("^\\s*cc\\s+", tnt_trailer, ignore.case = TRUE))) {
    warning("Columns (characters) were removed. Character indices in the trailing 'cc' (ordering) or NEXUS ASSUMPTIONS block are not automatically updated and may now be incorrect. Please review the output file.")
  }

  # Validate output_format
  if (!is.null(output_format)) {
    output_format <- tolower(output_format)
  }

  if (!is.null(output_format) && output_format == "tnt") {
    name <- paste0(output, "_FILTERED.tnt")

    # Enforce bracket format [01] for TNT
    data_tnt <- data
    data_tnt[] <- gsub("\\{|\\(", "[", data_tnt)
    data_tnt[] <- gsub("\\}|\\)", "]", data_tnt)

    if (is.null(tnt_trailer)) tnt_trailer <- "proc /;"

    out_lines <- character()
    if (length(tnt_preamble) > 0) out_lines <- c(out_lines, tnt_preamble)

    out_lines <- c(out_lines, "xread")
    out_lines <- c(out_lines, paste(ncol(data_tnt), nrow(data_tnt)))

    if (length(tnt_block_tag) > 0 && nzchar(tnt_block_tag)) {
      out_lines <- c(out_lines, tnt_block_tag)
    }

    for (i in seq_len(nrow(data_tnt))) {
      seq_str <- paste(data_tnt[i, ], collapse = "")
      out_lines <- c(out_lines, paste0(rownames(data_tnt)[i], " ", seq_str))
    }

    out_lines <- c(out_lines, ";")
    if (length(tnt_trailer) > 0) out_lines <- c(out_lines, tnt_trailer)

    writeLines(out_lines, name)
  } else {
    # Default to nexus
    name <- paste0(output, "_FILTERED.nexus")

    # Enforce curly brace format {01} for Nexus
    data_nexus <- data
    data_nexus[] <- gsub("\\[|\\(", "{", data_nexus)
    data_nexus[] <- gsub("\\]|\\)", "}", data_nexus)

    # Translate TNT ordering blocks (cc + ...) to NEXUS ASSUMPTIONS block
    assump_block <- character(0)
    if (!is.null(tnt_trailer)) {
      cc_lines <- grep("^\\s*cc\\s+\\+", tnt_trailer, ignore.case = TRUE, value = TRUE)
      if (length(cc_lines) > 0) {
        all_tokens <- character()
        for (line in cc_lines) {
          line_body <- sub("^\\s*cc\\s+\\+\\s+", "", line, ignore.case = TRUE)
          line_body <- sub(";.*$", "", line_body) # Remove semicolon and beyond
          tokens <- strsplit(trimws(line_body), "\\s+")[[1]]
          all_tokens <- c(all_tokens, tokens)
        }

        if (length(all_tokens) > 0) {
          # Establish old column index against newly dropped/surviving column index rules
          new_idx <- rep(NA, length(cols_to_delete))
          new_idx[!cols_to_delete] <- seq_len(sum(!cols_to_delete))

          # Dynamically expand and shift TNT 0-based indices to 1-based original index chunks
          idx_list <- lapply(all_tokens, function(x) {
            # Matches against classic TNT `.` interval format
            if (grepl("\\.", x)) {
              parts <- as.numeric(strsplit(x, "\\.+")[[1]])
              if (length(parts) == 2 && !any(is.na(parts))) {
                return(seq(parts[1] + 1, parts[2] + 1))
              }
            } else if (grepl("-", x)) {
              parts <- as.numeric(strsplit(x, "-")[[1]])
              if (length(parts) == 2 && !any(is.na(parts))) {
                return(seq(parts[1] + 1, parts[2] + 1))
              }
            }
            # Unspooled individual indices
            val <- as.numeric(x)
            if (!is.na(val)) {
              return(val + 1)
            }
            return(NULL)
          })

          # Remove non-assigned integers and flatten to vector array
          ordered_idx <- unlist(idx_list)

          if (length(ordered_idx) > 0) {
            # Shift characters depending on prior column removals, skipping over NAs effectively
            surviving_idx <- new_idx[ordered_idx]
            surviving_idx <- surviving_idx[!is.na(surviving_idx)]

            if (length(surviving_idx) > 0) {
              ordered_str <- paste(surviving_idx, collapse = " ")
              assump_block <- c(
                "BEGIN ASSUMPTIONS;",
                "    OPTIONS DEFTYPE=unord;",
                paste0("    TYPESET * Ordered = ", ordered_str, ";"),
                "END;"
              )
            }
          }
        }
      }
    }

    ape::write.nexus.data(data_nexus, file = name, format = "standard", interleaved = FALSE)
    temp <- readLines(name)
    temp <- gsub("INTERLEAVE=NO", "", temp)
    temp <- gsub("write.nexus.data.R", "RNODE", temp) # Clean Ape comment

    # Append the assumption block if one was constructed
    if (length(assump_block) > 0) {
      temp <- c(temp, assump_block)
    }

    writeLines(temp, name)
  }

  invisible(data)
}
