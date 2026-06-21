#' @title filterSharedTaxa
#' @name filterSharedTaxa
#' @description Filter two molecular or morphological matrices by shared and unique taxa (rows). Can output concatenated or separated matrices in NEXUS or TNT format.
#' @author Daniel YM Nakamura
#' @param input1 First input file (molecular or morphological matrix in nexus, fasta, or tnt format) or already loaded as a matrix object in R.
#' @param input1_format To load from a local file: 'nexus', 'fasta', or 'tnt'.
#' @param input2 Second input file (molecular or morphological matrix in nexus, fasta, or tnt format) or already loaded as a matrix object in R.
#' @param input2_format To load from a local file: 'nexus', 'fasta', or 'tnt'.
#' @param output_path Output path. For NEXUS output, "Desktop/Index" will output "Desktop/Index_SHARED.nexus" (concatenated) or "Desktop/Index_SHARED_1.nexus" and "Desktop/Index_SHARED_2.nexus" (separated). For TNT output, "Desktop/Index" will output "Desktop/Index.tnt" (concatenated) or "Desktop/Index_SHARED_1.tnt" and "Desktop/Index_SHARED_2.tnt" (separated). If NULL, returns matrices as a list without writing files.
#' @param output_format Output file format: 'nexus' or 'tnt'. If NULL (default), the format is inferred from output_path: paths ending in .tnt write TNT; all others write NEXUS.
#' @param output_concatenate Logical. If TRUE (default), concatenates the filtered matrices. If FALSE, returns them separately.
#' @param return_as_list Logical. If TRUE, returns the result as a list of matrices instead of writing to file. Default is FALSE.
#' @param level Taxon filtering level. Use 'strictly_shared' to keep only taxa shared by input1 and input2, 'shared+unique1' to keep shared taxa plus taxa unique to input1, or 'shared+unique2' to keep shared taxa plus taxa unique to input2. Missing data for taxa absent from one matrix are filled with '?'.
#' @examples
#' # Example with concatenated NEXUS output
#' filterSharedTaxa(input1="testdata/file1.nexus", input1_format="nexus",
#'                  input2="testdata/file2.nexus", input2_format="nexus",
#'                  output_path="testdata/shared", output_concatenate=TRUE)
#'
#' # Example with concatenated TNT output (protein data)
#' filterSharedTaxa(input1="testdata/013_MORPH_data.tnt", input1_format="tnt",
#'                  input2="testdata/013_MOL_data.tnt", input2_format="tnt",
#'                  output_path="testdata/013_TE_data", output_format="tnt",
#'                  output_concatenate=TRUE,
#'                  level="strictly_shared")
#'
#' # Example with separated output
#' filterSharedTaxa(input1="testdata/file1.nexus", input1_format="nexus",
#'                  input2="testdata/file2.nexus", input2_format="nexus",
#'                  output_path="testdata/shared", output_concatenate=FALSE)
#'
#' @export
filterSharedTaxa <- function(input1, input1_format,
                             input2, input2_format,
                             output_path = NULL,
                             output_format = NULL,
                             output_concatenate = TRUE,
                             return_as_list = FALSE,
                             level = "strictly_shared") {

  level <- match.arg(level, choices = c("strictly_shared", "shared+unique1", "shared+unique2"))
  if (!is.null(output_format)) {
    output_format <- match.arg(output_format, choices = c("nexus", "tnt"))
  }

  # Helper function to convert DNAbin to matrix with padding for unequal lengths
  dnabin_to_matrix <- function(dnabin_obj) {
    # Get sequence names
    names_seq <- names(dnabin_obj)

    # Find max length
    max_len <- max(sapply(dnabin_obj, length))

    # Convert each sequence to character and pad with 'N'
    mat <- sapply(dnabin_obj, function(x) {
      x_char <- as.character(x)
      if (length(x_char) < max_len) {
        x_char <- c(x_char, rep("N", max_len - length(x_char)))
      }
      x_char
    })

    # Transpose to get taxa as rows
    mat <- t(mat)
    rownames(mat) <- names_seq
    return(mat)
  }

  # Helper function to detect data type from nexus file
  detect_format_from_nexus <- function(nexus_file) {
    lines <- tolower(readLines(nexus_file))
    # Look for format declaration
    format_line <- grep("\\bformat\\b", lines, value = TRUE)
    if (length(format_line) > 0) {
      if (grepl("datatype\\s*=\\s*protein", format_line[1]) || grepl("protein", format_line[1])) {
        return("protein")
      } else if (grepl("datatype\\s*=\\s*dna", format_line[1]) || grepl("dna", format_line[1])) {
        return("dna")
      } else if (grepl("standard", format_line[1])) {
        return("standard")
      }
    }
    return("standard")
  }

  # Helper function to detect data type from TNT file
  detect_format_from_tnt <- function(tnt_file) {
    lines <- tolower(readLines(tnt_file, n = 10))  # Read first 10 lines
    # Look for protein state declarations
    if (any(grepl("nstates\\s+(prot|32)", lines))) {
      return("protein")
    }
    # Check if lines contain protein-specific amino acids in the data
    if (any(grepl("xread", lines))) {
      # If xread found, it's TNT format (standard is default for TNT)
      return("standard")
    }
    return("standard")
  }

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

    # Fast path for simple aligned sequences: protein/DNA/numeric matrices
    # without TNT polymorphism tokens.
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
        nrow = ntaxa,
        ncol = parsed_nchars,
        byrow = TRUE
      )
      rownames(mat) <- taxa_names
      return(mat)
    }

    # Slower, but still vectorized, path for polymorphic/inapplicable tokens.
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
    if (parsed_nchars != nchars) {
      warning("TNT header declares ", nchars, " characters, but parsed ",
              parsed_nchars, " characters in file: ", tnt_file)
    }

    mat <- matrix(unlist(mat_list, use.names = FALSE),
                  nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
    rownames(mat) <- taxa_names
    return(mat)
  }

  # Helper function to detect if matrix contains protein sequences
  detect_protein_from_matrix <- function(mat) {
    # Protein-specific amino acids (not in DNA or RNA)
    protein_only <- c("E", "F", "I", "K", "L", "M", "P", "Q", "Z")
    
    # Check if any protein-only amino acids are present
    mat_upper <- toupper(as.vector(mat))
    has_protein <- any(mat_upper %in% protein_only)
    
    # Also check for standard DNA ambiguity codes
    dna_codes <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-", "?")
    non_dna <- setdiff(unique(mat_upper[mat_upper != "-" & mat_upper != "?"]), c(dna_codes, protein_only))
    
    has_protein || (length(non_dna) > 0)
  }

  # Helper function to add rows of missing data for taxa absent from a matrix
  align_matrix_to_taxa <- function(mat, taxa) {
    aligned <- matrix("?", nrow = length(taxa), ncol = ncol(mat),
                      dimnames = list(taxa, colnames(mat)))
    present_taxa <- intersect(taxa, rownames(mat))
    aligned[present_taxa, ] <- mat[present_taxa, , drop = FALSE]
    aligned
  }

  # Helper function to catch invalid output directories before writing
  validate_output_files <- function(filenames) {
    output_dirs <- unique(dirname(filenames))
    output_dirs <- output_dirs[output_dirs != "."]

    missing_dirs <- output_dirs[!dir.exists(output_dirs)]
    if (length(missing_dirs) > 0) {
      msg <- paste0("Output directory does not exist: ",
                    paste(missing_dirs, collapse = ", "))
      if (any(grepl("^/", filenames))) {
        msg <- paste0(msg,
                      "\nPaths starting with '/' are absolute paths from the filesystem root. ",
                      "Use 'tnt/056_', './tnt/056_', or the full project path if you mean a local folder.")
      }
      stop(msg, call. = FALSE)
    }

    unwritable_dirs <- output_dirs[file.access(output_dirs, 2) != 0]
    if (length(unwritable_dirs) > 0) {
      stop("Output directory is not writable: ",
           paste(unwritable_dirs, collapse = ", "), call. = FALSE)
    }
  }

  # Helper function to write TNT format
  write_tnt_lines <- function(lines, filename) {
    if (length(lines) > 0 && lines[1] == "nstates 32") {
      text <- paste0(lines[1], "\n", paste(lines[-1], collapse = "\r\n"), "\r\n")
      writeBin(charToRaw(text), filename)
    } else {
      writeLines(lines, filename, sep = "\r\n")
    }
  }

  # Helper function to write TNT format
  write_tnt_format <- function(mat, filename, is_protein = FALSE) {
    ntaxa <- nrow(mat)
    nchars <- ncol(mat)
    
    # Create TNT output lines (header is nchars ntaxa, NOT ntaxa nchars)
    lines <- c("xread", paste(nchars, ntaxa))
    if (is_protein) {
      lines <- c("nstates 32", lines)
    }
    
    # Add protein block header if needed
    if (is_protein) {
      lines <- c(lines, "", "& [prot]")
    } else {
      lines <- c(lines, "", "& [num]")
    }
    
    # Add sequences
    for (i in 1:ntaxa) {
      taxon <- rownames(mat)[i]
      sequence <- paste(mat[i, ], collapse = "")
      lines <- c(lines, paste0(taxon, "\t", sequence))
    }
    
    # Add footer
    lines <- c(lines, ";", "", "", "proc /;", "comments 0", ";", "", "")
    
    # Write to file with TNT line endings
    write_tnt_lines(lines, filename)
  }

  # Helper function to write TNT format with multiple blocks (for concatenation)
  write_tnt_format_concat <- function(mat1, mat2, filename, format1 = "standard", format2 = "standard") {
    ntaxa <- nrow(mat1)
    nchars1 <- ncol(mat1)
    nchars2 <- ncol(mat2)
    total_nchars <- nchars1 + nchars2
    
    # Create TNT output lines (header is nchars ntaxa, NOT ntaxa nchars)
    lines <- c("xread", paste(total_nchars, ntaxa))
    if (format1 == "protein" || format2 == "protein") {
      lines <- c("nstates 32", lines)
    }
    
    # Determine block headers
    block1_header <- if (format1 == "protein") "& [prot]" else "& [num]"
    block2_header <- if (format2 == "protein") "& [prot]" else "& [num]"
    
    # Reorder blocks: protein should come first, then standard
    if (format1 == "protein" && format2 != "protein") {
      # mat1 is protein, mat2 is not - write in correct order
      lines <- c(lines, "", block1_header)
      for (i in 1:ntaxa) {
        taxon <- rownames(mat1)[i]
        sequence <- paste(mat1[i, ], collapse = "")
        lines <- c(lines, paste0(taxon, "\t", sequence))
      }
      lines <- c(lines, "", block2_header)
      for (i in 1:ntaxa) {
        taxon <- rownames(mat2)[i]
        sequence <- paste(mat2[i, ], collapse = "")
        lines <- c(lines, paste0(taxon, "\t", sequence))
      }
    } else if (format1 != "protein" && format2 == "protein") {
      # mat1 is not protein, mat2 is protein - swap order for output
      lines <- c(lines, "", block2_header)
      for (i in 1:ntaxa) {
        taxon <- rownames(mat2)[i]
        sequence <- paste(mat2[i, ], collapse = "")
        lines <- c(lines, paste0(taxon, "\t", sequence))
      }
      lines <- c(lines, "", block1_header)
      for (i in 1:ntaxa) {
        taxon <- rownames(mat1)[i]
        sequence <- paste(mat1[i, ], collapse = "")
        lines <- c(lines, paste0(taxon, "\t", sequence))
      }
    } else {
      # Both protein or both standard - write in input order
      lines <- c(lines, "", block1_header)
      for (i in 1:ntaxa) {
        taxon <- rownames(mat1)[i]
        sequence <- paste(mat1[i, ], collapse = "")
        lines <- c(lines, paste0(taxon, "\t", sequence))
      }
      lines <- c(lines, "", block2_header)
      for (i in 1:ntaxa) {
        taxon <- rownames(mat2)[i]
        sequence <- paste(mat2[i, ], collapse = "")
        lines <- c(lines, paste0(taxon, "\t", sequence))
      }
    }
    
    # Add footer
    lines <- c(lines, ";", "", "", "proc /;", "comments 0", ";", "", "")
    
    # Write to file with TNT line endings
    write_tnt_lines(lines, filename)
  }

  # Track format types for concatenation
  format1 <- "standard"
  format2 <- "standard"

  # Load first matrix
  if (is.matrix(input1) || is.data.frame(input1)) {
    mat1 <- as.matrix(input1)
  } else {
    if (is.null(input1_format)) {
      stop("If 'input1' is a file path, you must provide input1_format = 'nexus', 'fasta', or 'tnt'")
    }
    if (input1_format == "nexus") {
      format1 <- detect_format_from_nexus(input1)
      mat1 <- TreeTools::ReadCharacters(input1)
    } else if (input1_format == "fasta") {
      mat1 <- ape::read.FASTA(input1)
      # Handle DNAbin objects with unequal sequence lengths
      if (inherits(mat1, "DNAbin")) {
        mat1 <- dnabin_to_matrix(mat1)
        format1 <- "dna"
      } else {
        mat1 <- as.matrix(mat1)
        # Try to detect if it's protein
        if (detect_protein_from_matrix(mat1)) {
          format1 <- "protein"
        }
      }
    } else if (input1_format == "tnt") {
      format1 <- detect_format_from_tnt(input1)
      mat1 <- read_tnt_manual(input1)
    } else {
      stop("input1_format must be 'nexus', 'fasta', or 'tnt'")
    }
    mat1 <- as.matrix(mat1)
  }

  # Load second matrix
  if (is.matrix(input2) || is.data.frame(input2)) {
    mat2 <- as.matrix(input2)
  } else {
    if (is.null(input2_format)) {
      stop("If 'input2' is a file path, you must provide input2_format = 'nexus', 'fasta', or 'tnt'")
    }
    if (input2_format == "nexus") {
      format2 <- detect_format_from_nexus(input2)
      mat2 <- TreeTools::ReadCharacters(input2)
    } else if (input2_format == "fasta") {
      mat2 <- ape::read.FASTA(input2)
      # Handle DNAbin objects with unequal sequence lengths
      if (inherits(mat2, "DNAbin")) {
        mat2 <- dnabin_to_matrix(mat2)
        format2 <- "dna"
      } else {
        mat2 <- as.matrix(mat2)
        # Try to detect if it's protein
        if (detect_protein_from_matrix(mat2)) {
          format2 <- "protein"
        }
      }
    } else if (input2_format == "tnt") {
      format2 <- detect_format_from_tnt(input2)
      mat2 <- read_tnt_manual(input2)
    } else {
      stop("input2_format must be 'nexus', 'fasta', or 'tnt'")
    }
    mat2 <- as.matrix(mat2)
  }

  # Get taxa names
  taxa1 <- rownames(mat1)
  taxa2 <- rownames(mat2)

  if (is.null(taxa1) || is.null(taxa2)) {
    stop("Both matrices must have row names (taxa names)")
  }

  # Find shared taxa
  shared_taxa <- intersect(taxa1, taxa2)
  unique_taxa1 <- setdiff(taxa1, taxa2)
  unique_taxa2 <- setdiff(taxa2, taxa1)

  selected_taxa <- switch(level,
                          strictly_shared = shared_taxa,
                          `shared+unique1` = taxa1,
                          `shared+unique2` = taxa2)

  if (length(selected_taxa) == 0) {
    stop("No taxa selected with the requested level")
  }

  cat("Total taxa in matrix 1:", length(taxa1), "\n")
  cat("Total taxa in matrix 2:", length(taxa2), "\n")
  cat("Shared taxa:", length(shared_taxa), "\n")
  cat("Unique taxa in matrix 1:", length(unique_taxa1), "\n")
  cat("Unique taxa in matrix 2:", length(unique_taxa2), "\n")
  cat("Level:", level, "\n")
  cat("Selected taxa:", length(selected_taxa), "\n")
  cat("Shared taxa:", paste(shared_taxa, collapse = ", "), "\n\n")

  # Filter matrices to keep the selected taxa, filling absent rows with missing data
  mat1_filtered <- align_matrix_to_taxa(mat1, selected_taxa)
  mat2_filtered <- align_matrix_to_taxa(mat2, selected_taxa)

  # Return as list if requested
  if (return_as_list) {
    if (output_concatenate) {
      result_mat <- cbind(mat1_filtered, mat2_filtered)
      return(list(concatenated = result_mat))
    } else {
      return(list(matrix1 = mat1_filtered, matrix2 = mat2_filtered))
    }
  }

  # Write to file if output_path is provided
  if (!is.null(output_path)) {
    # Determine if output should be TNT format
    if (is.null(output_format)) {
      output_file_format <- if (grepl("\\.tnt$", output_path, ignore.case = TRUE)) "tnt" else "nexus"
    } else {
      output_file_format <- output_format
    }
    is_tnt_output <- output_file_format == "tnt"
    
    if (output_concatenate) {
      # Concatenate matrices
      result_mat <- cbind(mat1_filtered, mat2_filtered)
      
      if (is_tnt_output) {
        output_name <- output_path
        if (!grepl("\\.tnt$", output_name, ignore.case = TRUE)) {
          output_name <- paste0(output_path, ".tnt")
        }
      } else {
        output_name <- paste0(output_path, "_SHARED.nexus")
      }
      validate_output_files(output_name)

      # Determine NEXUS data type (use protein if either input is protein, DNA if either is DNA, else standard)
      nexus_data_format <- "standard"
      if (format1 == "protein" || format2 == "protein") {
        nexus_data_format <- "protein"
      } else if (format1 == "dna" || format2 == "dna") {
        nexus_data_format <- "dna"
      }

      if (is_tnt_output) {
        # Write TNT format with separate blocks for each matrix
        write_tnt_format_concat(mat1_filtered, mat2_filtered, output_name, format1 = format1, format2 = format2)
        cat("Concatenated matrix written to:", output_name, "\n")
        cat("Format: TNT (blocks: ", format1, " +", format2, ")\n")
      } else {
        # Write NEXUS format
        ape::write.nexus.data(result_mat, file = output_name, format = nexus_data_format, interleaved = F)
        temp <- gsub("INTERLEAVE=NO", "", readLines(output_name))
        temp <- gsub("write.nexus.data.R", "RNODE", temp)
        writeLines(temp, output_name)
        cat("Concatenated matrix written to:", output_name, "\n")
        cat("Format:", nexus_data_format, "\n")
      }
    } else {
      # Separate outputs
      if (is_tnt_output) {
        # Extract base path for TNT files
        base_path <- sub("\\.tnt$", "", output_path, ignore.case = TRUE)
        output_name1 <- paste0(base_path, "_SHARED_1.tnt")
        output_name2 <- paste0(base_path, "_SHARED_2.tnt")
        validate_output_files(c(output_name1, output_name2))
        
        is_protein1 <- (format1 == "protein")
        is_protein2 <- (format2 == "protein")
        
        write_tnt_format(mat1_filtered, output_name1, is_protein = is_protein1)
        write_tnt_format(mat2_filtered, output_name2, is_protein = is_protein2)
      } else {
        output_name1 <- paste0(output_path, "_SHARED_1.nexus")
        output_name2 <- paste0(output_path, "_SHARED_2.nexus")
        validate_output_files(c(output_name1, output_name2))

        ape::write.nexus.data(mat1_filtered, file = output_name1, format = format1, interleaved = F)
        temp1 <- gsub("INTERLEAVE=NO", "", readLines(output_name1))
        temp1 <- gsub("write.nexus.data.R", "RNODE", temp1)
        writeLines(temp1, output_name1)

        ape::write.nexus.data(mat2_filtered, file = output_name2, format = format2, interleaved = F)
        temp2 <- gsub("INTERLEAVE=NO", "", readLines(output_name2))
        temp2 <- gsub("write.nexus.data.R", "RNODE", temp2)
        writeLines(temp2, output_name2)
      }

      cat("Filtered matrices written to:\n")
      cat("  -", output_name1, "(format:", format1, ")\n")
      cat("  -", output_name2, "(format:", format2, ")\n")
    }
  }

  # Return invisibly the filtered matrices
  if (output_concatenate) {
    result_mat <- cbind(mat1_filtered, mat2_filtered)
    return(invisible(list(concatenated = result_mat)))
  } else {
    return(invisible(list(matrix1 = mat1_filtered, matrix2 = mat2_filtered)))
  }
}
