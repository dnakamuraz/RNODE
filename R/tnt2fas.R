#' @title tnt2fas
#' @name tnt2fas
#' @description Converts a TNT matrix into FASTA format. Only works for DNA and PROTEIN data.
#' @author Daniel YM Nakamura
#' @param input_file Path to the input TNT file.
#' @param output_file Path to the output FASTA file.
#' @param iupac_ambiguity Logical; if TRUE, uses IUPAC symbols for ambiguities (default TRUE for FASTA).
#' @export
tnt2fas <- function(input_file, output_file, iupac_ambiguity = TRUE) {
  lines <- readLines(input_file, warn = FALSE)
  datatype <- "STANDARD"
  in_matrix <- FALSE
  taxa_seqs <- list()

  parse_xread_dims <- function(line) {
    cleaned <- gsub("['\"]", "", line)
    nums <- unlist(regmatches(cleaned, gregexpr("\\b[0-9]+\\b", cleaned, perl = TRUE)))
    if (length(nums) >= 2) {
      return(as.numeric(nums[1:2]))
    }
    return(NULL)
  }

  i <- 1
  while (i <= length(lines)) {
    line <- trimws(lines[i])
    line_lower <- tolower(line)

    if (line_lower == "xread" || grepl("^xread", line_lower) || grepl("^nstates", line_lower)) {
      if (line_lower == "xread" && i < length(lines)) {
        i <- i + 1
        line <- trimws(lines[i])
      }
      dims <- parse_xread_dims(line)
      if (is.null(dims)) {
        i <- i + 1
        if (i <= length(lines)) {
          dims <- parse_xread_dims(lines[i])
        }
      }
      in_matrix <- TRUE
      i <- i + 1
      next
    }

    if (in_matrix) {
      if (line == ";" || line_lower == "proc /;") {
        in_matrix <- FALSE
      } else if (startsWith(line_lower, "&[dna") || startsWith(line_lower, "& [dna")) {
        datatype <- "DNA"
      } else if (startsWith(line_lower, "&[prot") || startsWith(line_lower, "& [prot")) {
        datatype <- "PROTEIN"
      } else if (startsWith(line_lower, "&[num") || startsWith(line_lower, "& [num")) {
        datatype <- "STANDARD"
      } else if (line != "") {
        if (datatype %in% c("DNA", "PROTEIN")) {
          parts <- unlist(strsplit(line, "\\s+"))
          if (length(parts) >= 2) {
            taxon <- parts[1]
            seq_data <- paste(parts[-1], collapse = "")
            
            if (iupac_ambiguity) {
              matches <- gregexpr("\\[([^\\]]+)\\]", seq_data)
              ambigs <- regmatches(seq_data, matches)[[1]]
              if (length(ambigs) > 0) {
                for (ambig in unique(ambigs)) {
                  chars <- gsub("[\\[\\]]", "", ambig)
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
                  
                  pattern <- gsub("([\\[\\]])", "\\\\\\1", ambig)
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
    i <- i + 1
  }

  if (length(taxa_seqs) == 0) {
    stop("No DNA or PROTEIN sequences parsed. Check input file format.")
  }

  out_lines <- c()
  for (taxon in names(taxa_seqs)) {
    out_lines <- c(out_lines, paste0(">", taxon), taxa_seqs[[taxon]])
  }
  
  writeLines(out_lines, output_file)
  message(paste("Successfully converted to", output_file))
}
