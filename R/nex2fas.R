#' @title nex2fas
#' @name nex2fas
#' @description Converts a NEXUS matrix into FASTA format. Only works for DNA and PROTEIN data.
#' @author Daniel YM Nakamura
#' @param input_file Path to the input NEXUS file.
#' @param output_file Path to the output FASTA file.
#' @param iupac_ambiguity Logical; if TRUE, uses IUPAC symbols for ambiguities (default TRUE for FASTA).
#' @export
nex2fas <- function(input_file, output_file, iupac_ambiguity = TRUE) {
  lines <- readLines(input_file, warn = FALSE)
  datatype <- "STANDARD"
  in_matrix <- FALSE
  taxa_seqs <- list()

  for (i in seq_along(lines)) {
    line <- trimws(lines[i])
    line_upper <- toupper(line)

    if (grepl("FORMAT", line_upper) && grepl("DATATYPE", line_upper)) {
      if (grepl("DATATYPE\\s*=\\s*DNA", line_upper)) {
        datatype <- "DNA"
      } else if (grepl("DATATYPE\\s*=\\s*PROTEIN", line_upper)) {
        datatype <- "PROTEIN"
      } else if (grepl("DATATYPE\\s*=\\s*STANDARD", line_upper)) {
        datatype <- "STANDARD"
      }
    }

    if (line_upper == "MATRIX") {
      in_matrix <- TRUE
      if (datatype == "STANDARD") {
        stop("nex2fas only works for DNA and PROTEIN data.")
      }
      next
    }

    if (in_matrix) {
      if (line_upper == ";" || grepl("^END;", line_upper)) {
        in_matrix <- FALSE
      } else if (line != "") {
        if (!startsWith(line, "[")) {
          parts <- unlist(strsplit(line, "\\s+"))
          if (length(parts) >= 2) {
            taxon <- parts[1]
            seq_data <- paste(parts[-1], collapse = "")
            
            if (iupac_ambiguity) {
              matches <- gregexpr("[\\{\\(]([^\\}\\)]+)[\\}\\)]", seq_data)
              ambigs <- regmatches(seq_data, matches)[[1]]
              if (length(ambigs) > 0) {
                for (ambig in unique(ambigs)) {
                  chars <- gsub("[\\(\\)\\{\\}]", "", ambig)
                  chars <- toupper(strsplit(chars, "")[[1]])
                  chars <- chars[chars != " " & chars != ""]
                  chars <- sort(unique(chars))
                  chars_str <- paste(chars, collapse = "")
                  
                  iupac_code <- switch(chars_str,
                                       "AG" = "R", "CT" = "Y", "CG" = "S", "AT" = "W",
                                       "GT" = "K", "AC" = "M", "CGT" = "B", "AGT" = "D",
                                       "ACT" = "H", "ACG" = "V", "ACGT" = "N",
                                       "DN" = "B", "EQ" = "Z", "IL" = "J",
                                       "ACDEFGHIKLMNPQRSTVWY" = "X", ambig)
                  
                  pattern <- gsub("([\\(\\)\\{\\}])", "\\\\\\1", ambig)
                  seq_data <- gsub(pattern, iupac_code, seq_data)
                }
              }
            }
            
            if (is.null(taxa_seqs[[taxon]])) {
              taxa_seqs[[taxon]] <- seq_data
            } else {
              taxa_seqs[[taxon]] <- paste0(taxa_seqs[[taxon]], seq_data)
            }
          }
        }
      }
    }
  }

  if (length(taxa_seqs) == 0) {
    stop("No sequences parsed. Check input file format.")
  }

  out_lines <- c()
  for (taxon in names(taxa_seqs)) {
    out_lines <- c(out_lines, paste0(">", taxon), taxa_seqs[[taxon]])
  }
  
  writeLines(out_lines, output_file)
  message(paste("Successfully converted to", output_file))
}
