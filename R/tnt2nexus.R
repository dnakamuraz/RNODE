#' @title tnt2nexus
#' @name tnt2nexus
#' @description Converts a TNT matrix into NEXUS format while preserving header data type and ordered character structures.
#' @author Daniel YM Nakamura
#' @param input_file Path to the input TNT file.
#' @param output_file Path to the output NEXUS file.
#' @export
tnt2nexus <- function(input_file, output_file) {
  lines <- readLines(input_file, warn = FALSE)

  nchar <- 0
  ntax <- 0
  datatype <- "STANDARD"
  matrix_lines <- c()
  ordered_chars <- c()

  in_matrix <- FALSE

  i <- 1
  while (i <= length(lines)) {
    line <- trimws(lines[i])
    line_lower <- tolower(line)

    if (line == "" || startsWith(line, "'")) {
      i <- i + 1
      next
    }

    if (line_lower == "xread") {
      # read next non-empty line for dims
      i <- i + 1
      while (i <= length(lines) && (trimws(lines[i]) == "" || startsWith(trimws(lines[i]), "'"))) {
        i <- i + 1
      }
      if (i <= length(lines)) {
        dims <- unlist(strsplit(trimws(lines[i]), "\\s+"))
        if (length(dims) >= 2) {
          nchar <- as.numeric(dims[1])
          ntax <- as.numeric(dims[2])
        }
      }

      # start matrix read
      in_matrix <- TRUE
      i <- i + 1
      next
    }

    if (in_matrix) {
      if (line == ";" || line_lower == "proc /;") {
        in_matrix <- FALSE
      } else if (startsWith(line_lower, "&[dna]") || startsWith(line_lower, "& [dna]")) {
        datatype <- "DNA"
      } else if (startsWith(line_lower, "&[prot]") || startsWith(line_lower, "& [prot]")) {
        datatype <- "PROTEIN"
      } else if (startsWith(line_lower, "&[num]") || startsWith(line_lower, "& [num]")) {
        datatype <- "STANDARD"
      } else if (line != "") {
        matrix_lines <- c(matrix_lines, line)
      }
      i <- i + 1
      next
    }

    # Check for ccode +
    if (startsWith(line_lower, "ccode +")) {
      ord_part <- gsub("ccode \\+\\s*", "", line_lower)
      ord_part <- gsub(";", "", ord_part)

      tokens <- unlist(strsplit(trimws(ord_part), "\\s+"))
      for (tok in tokens) {
        if (grepl("\\.", tok)) {
          pts <- as.numeric(unlist(strsplit(tok, "\\.")))
          if (!is.na(pts[1]) && !is.na(pts[2])) {
            ordered_chars <- c(ordered_chars, pts[1]:pts[2])
          }
        } else if (grepl("-", tok)) {
          pts <- as.numeric(unlist(strsplit(tok, "-")))
          if (!is.na(pts[1]) && !is.na(pts[2])) {
            ordered_chars <- c(ordered_chars, pts[1]:pts[2])
          }
        } else {
          if (tok != "") {
            val <- as.numeric(tok)
            if (!is.na(val)) ordered_chars <- c(ordered_chars, val)
          }
        }
      }
    }

    i <- i + 1
  }

  # Format NEXUS output
  out_lines <- c(
    "#NEXUS", "BEGIN DATA;",
    paste0("\tDIMENSIONS NTAX=", ntax, " NCHAR=", nchar, ";")
  )

  if (datatype == "STANDARD") {
    out_lines <- c(out_lines, paste0("\tFORMAT DATATYPE=", datatype, " MISSING=? GAP=- SYMBOLS=\"0 1 2 3 4 5 6 7 8 9\";"))
  } else {
    out_lines <- c(out_lines, paste0("\tFORMAT DATATYPE=", datatype, " MISSING=? GAP=-;"))
  }

  out_lines <- c(out_lines, "\tMATRIX")
  for (ml in matrix_lines) {
    out_lines <- c(out_lines, paste0("\t\t", ml))
  }
  out_lines <- c(out_lines, "\t;", "END;", "")

  if (length(ordered_chars) > 0) {
    # Convert 0-based TNT to 1-based NEXUS
    ordered_chars <- unique(ordered_chars) + 1
    ordered_str <- paste(ordered_chars, collapse = " ")
    out_lines <- c(
      out_lines, "BEGIN ASSUMPTIONS;",
      paste0("\tTYPESET * UNTITLED = ord: ", ordered_str, ";"),
      "END;"
    )
  }

  writeLines(out_lines, output_file)
  message(paste("Successfully converted to", output_file))
}
