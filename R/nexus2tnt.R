#' @title nexus2tnt
#' @name nexus2tnt
#' @description Converts a NEXUS matrix into TNT format while preserving header data type and ordered character structures.
#' @author Daniel YM Nakamura
#' @param input_file Path to the input NEXUS file.
#' @param output_file Path to the output TNT file.
#' @export
nexus2tnt <- function(input_file, output_file) {
  lines <- readLines(input_file, warn = FALSE)

  ntax <- 0
  nchar <- 0
  datatype <- "STANDARD"
  matrix_lines <- c()
  ordered_chars <- c()

  in_matrix <- FALSE
  in_assumptions <- FALSE

  for (i in seq_along(lines)) {
    line <- trimws(lines[i])
    line_upper <- toupper(line)

    # Dimensions
    if (grepl("NTAX\\s*=", line_upper)) {
      ntax <- as.integer(gsub(".*NTAX\\s*=\\s*([0-9]+).*", "\\1", line_upper))
    }
    if (grepl("NCHAR\\s*=", line_upper)) {
      nchar <- as.integer(gsub(".*NCHAR\\s*=\\s*([0-9]+).*", "\\1", line_upper))
    }

    # Format
    if (grepl("FORMAT", line_upper) && grepl("DATATYPE", line_upper)) {
      if (grepl("DATATYPE\\s*=\\s*DNA", line_upper)) {
        datatype <- "DNA"
      } else if (grepl("DATATYPE\\s*=\\s*PROTEIN", line_upper)) {
        datatype <- "PROTEIN"
      } else if (grepl("DATATYPE\\s*=\\s*STANDARD", line_upper)) {
        datatype <- "STANDARD"
      }
    }

    # Matrix parsing
    if (line_upper == "MATRIX") {
      in_matrix <- TRUE
      next
    }

    if (in_matrix) {
      if (line_upper == ";" || grepl("^END;", line_upper)) {
        in_matrix <- FALSE
      } else if (line != "") {
        # Check if line contains taxon data (starts with alphanumeric or quote, not [)
        if (!startsWith(line, "[")) {
          parts <- unlist(strsplit(line, "\\s+"))
          if (length(parts) >= 2) {
            taxon <- parts[1]
            seq_data <- paste(parts[-1], collapse = "")
            seq_data <- gsub("\\(|\\{", "[", seq_data)
            seq_data <- gsub("\\)|\\}", "]", seq_data)
            matrix_lines <- c(matrix_lines, paste(taxon, seq_data))
          }
        }
      }
    }

    # Assumptions parsing (ordered characters)
    if (line_upper == "BEGIN ASSUMPTIONS;") {
      in_assumptions <- TRUE
      next
    }
    if (in_assumptions) {
      if (line_upper == "END;") {
        in_assumptions <- FALSE
      } else if (grepl("TYPESET.*=.*ORD:", line_upper)) {
        # Extract everything after "ord:"
        ord_part <- gsub(".*ORD:\\s*(.*)", "\\1", line_upper)
        ord_part <- gsub(";", "", ord_part)

        # Remove all spaces to handle format like "1 - 193"
        s_nospc <- gsub("\\s+", "", ord_part)
        # Parse ranges separated by comma
        parts <- unlist(strsplit(s_nospc, ","))
        for (p in parts) {
          if (grepl("-", p)) {
            b <- as.numeric(unlist(strsplit(p, "-")))
            if (!is.na(b[1]) && !is.na(b[2])) {
              ordered_chars <- c(ordered_chars, b[1]:b[2])
            }
          } else {
            if (p != "") {
              ordered_chars <- c(ordered_chars, as.numeric(p))
            }
          }
        }
      }
    }
  }

  # Format TNT output
  if (datatype == "PROTEIN") {
    tnt_header <- c("nstates 32", "xread")
  } else {
    tnt_header <- "xread"
  }

  tnt_datatype <- "&[num]"
  if (datatype == "DNA") {
    tnt_datatype <- "&[dna]"
  } else if (datatype == "PROTEIN") {
    tnt_datatype <- "&[prot]"
  }

  out_lines <- c(tnt_header, paste(nchar, ntax), "", tnt_datatype)
  out_lines <- c(out_lines, matrix_lines, ";")

  if (length(ordered_chars) > 0) {
    ordered_chars <- unique(ordered_chars)
    # Convert from 1-based (NEXUS) to 0-based (TNT)
    tnt_ordered <- ordered_chars - 1
    tnt_ordered_str <- paste(tnt_ordered, collapse = " ")
    out_lines <- c(out_lines, paste("ccode +", tnt_ordered_str, ";"))
  }

  out_lines <- c(out_lines, "", "proc /;", "comments 0", ";")

  writeLines(out_lines, output_file)
  message(paste("Successfully converted to", output_file))
}
