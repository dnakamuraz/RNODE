#' @title filterTaxaList
#' @name filterTaxaList
#' @description Filter a molecular or morphological matrix by keeping or deleting specific taxa based on a provided list.
#' @author Daniel YM Nakamura
#' @param input A local file path (molecular or morphological matrix in nexus, tnt, or fasta format) or already loaded as a matrix object in R.
#' @param input_format To load from a local file: 'nexus', 'tnt', 'fasta', or NULL. Defaults to NULL (auto-detects from extension).
#' @param list A local .txt file path or an R object (vector, list, or data.frame) containing the list of terminals, with one terminal name per row/element.
#' @param filter Character string. Use 'delete' to remove taxa contained in the list, or 'keep' to retain only the taxa contained in the list.
#' @param output Output path (e.g., "Desktop/Filtered"). If NULL, returns the matrix as an R object without writing files.
#' @param output_format Output file format: 'nexus', 'tnt', 'fasta', or NULL. If NULL, the format is inferred from the `output` extension, or simply returns as an R object if `output` is also NULL. FASTA output is only available for DNA and protein data.
#' @examples
#' \dontrun{
#' # Keep only the taxa specified in a text file
#' filterTaxaList(
#'     input = "testdata/matrix.tnt",
#'     list = "testdata/taxa_to_keep.txt",
#'     filter = "keep",
#'     output = "testdata/filtered_matrix.tnt"
#' )
#'
#' # Delete taxa using an R vector and return as an object
#' taxa_to_remove <- c("TaxonA", "TaxonB")
#' filtered_mat <- filterTaxaList(
#'     input = "testdata/matrix.nexus",
#'     list = taxa_to_remove,
#'     filter = "delete"
#' )
#'
#' # Filter a FASTA alignment and write the result back to FASTA
#' filterTaxaList(
#'     input = "testdata/alignment.fasta",
#'     input_format = "fasta",
#'     list = "testdata/taxa_to_keep.txt",
#'     filter = "keep",
#'     output = "testdata/filtered_alignment.fasta"
#' )
#' }
#' @export
filterTaxaList <- function(input, input_format = NULL, list, filter = "keep",
                           output = NULL, output_format = NULL) {
    filter <- match.arg(filter, choices = c("delete", "keep"))

    if (!is.null(input_format)) {
        input_format <- match.arg(input_format, choices = c("nexus", "tnt", "fasta"))
    }

    if (!is.null(output_format)) {
        output_format <- match.arg(output_format, choices = c("nexus", "tnt", "fasta"))
    }

    # Helper functions (Reused from RNODE standards)
    detect_format_from_nexus <- function(nexus_file) {
        lines <- tolower(readLines(nexus_file))
        format_line <- grep("\\bformat\\b", lines, value = TRUE)
        if (length(format_line) > 0) {
            if (grepl("datatype\\s*=\\s*protein", format_line[1]) || grepl("protein", format_line[1])) {
                return("protein")
            }
            if (grepl("datatype\\s*=\\s*dna", format_line[1]) || grepl("dna", format_line[1])) {
                return("dna")
            }
            if (grepl("datatype\\s*=\\s*rna", format_line[1]) || grepl("rna", format_line[1])) {
                return("dna")
            }
            if (grepl("datatype\\s*=\\s*standard", format_line[1]) || grepl("standard", format_line[1])) {
                return("standard")
            }
        }
        return("unknown")
    }

    detect_format_from_tnt <- function(tnt_file) {
        lines <- readLines(tnt_file)
        lines_lc <- tolower(lines)
        if (any(grepl("nstates\\s+(prot|32)", lines_lc))) {
            return("protein")
        }
        if (any(grepl("nstates\\s+dna", lines_lc))) {
            return("dna")
        }
        if (any(grepl("xread\\s+dna", lines_lc))) {
            return("dna")
        }
        if (any(grepl("^\\s*&\\s*\\[dna\\]", lines_lc))) {
            return("dna")
        }
        if (any(grepl("^\\s*&\\s*\\[prot\\]", lines_lc))) {
            return("protein")
        }
        if (any(grepl("^\\s*&\\s*\\[num\\]", lines_lc))) {
            return("standard")
        }
        return("unknown")
    }

    detect_format_from_matrix <- function(mat) {
        mat_upper <- toupper(as.vector(mat))
        mat_upper <- gsub("[][(){}]", "", mat_upper)
        chars <- unlist(strsplit(mat_upper, ""))
        valid_chars <- chars[!(chars %in% c("?", "-", "N", "X", ""))]
        if (length(valid_chars) == 0) {
            return("standard")
        }
        if (any(grepl("[0-9]", valid_chars))) {
            return("standard")
        }

        protein_only <- c("E", "F", "I", "L", "P", "Q", "Z", "J", "O")
        if (any(valid_chars %in% protein_only)) {
            return("protein")
        }

        dna_codes <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")
        if (all(valid_chars %in% dna_codes)) {
            return("dna")
        }
        if (all(grepl("^[A-Z]$", valid_chars))) {
            return("protein")
        }
        return("standard")
    }

    detect_format_from_fasta <- function(fasta_file) {
        lines <- readLines(fasta_file, warn = FALSE)
        seq_lines <- lines[!startsWith(trimws(lines), ">")]
        seq_lines <- gsub("\\s+", "", seq_lines)
        seq_chars <- unique(strsplit(toupper(paste(seq_lines, collapse = "")), "")[[1]])
        seq_chars <- seq_chars[nzchar(seq_chars)]
        if (length(seq_chars) == 0) {
            return("unknown")
        }

        dna_codes <- c(
            "A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M",
            "B", "D", "H", "V", "N", "-", "?", "."
        )
        if (all(seq_chars %in% dna_codes)) {
            return("dna")
        }

        protein_codes <- c(
            "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L",
            "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z",
            "X", "J", "U", "O", "*", "-", "?", "."
        )
        if (all(seq_chars %in% protein_codes)) {
            return("protein")
        }
        return("unknown")
    }

    fasta_to_matrix <- function(fasta_file, data_format = "dna") {
        fasta_type <- if (identical(data_format, "protein")) "AA" else "DNA"
        fasta_obj <- ape::read.FASTA(fasta_file, type = fasta_type)
        seq_list <- as.character(fasta_obj)

        if (!is.list(seq_list)) {
            seq_list <- list(seq_list)
            names(seq_list) <- names(fasta_obj)
        }

        pad <- if (identical(data_format, "protein")) "X" else "N"
        max_len <- max(lengths(seq_list))
        mat <- t(vapply(seq_list, function(seq_chars) {
            seq_chars <- toupper(seq_chars)
            if (length(seq_chars) < max_len) {
                seq_chars <- c(seq_chars, rep(pad, max_len - length(seq_chars)))
            }
            seq_chars
        }, character(max_len)))
        rownames(mat) <- names(seq_list)
        return(mat)
    }

    fix_ambiguity_for_fasta <- function(sequence, data_format = "dna") {
        tokens <- regmatches(
            sequence,
            gregexpr("[\\{\\(\\[]([^\\}\\)\\]]+)[\\}\\)\\]]", sequence, perl = TRUE)
        )[[1]]

        if (length(tokens) > 0) {
            for (token in unique(tokens)) {
                chars <- gsub("[^A-Za-z]", "", token)
                chars <- sort(unique(strsplit(toupper(chars), "")[[1]]))
                chars_str <- paste(chars, collapse = "")
                if (identical(data_format, "protein")) {
                    iupac <- switch(chars_str,
                        "DN" = "B", "EQ" = "Z", "IL" = "J",
                        "X"
                    )
                } else {
                    iupac <- switch(chars_str,
                        "AG" = "R", "CT" = "Y", "CG" = "S", "AT" = "W",
                        "GT" = "K", "AC" = "M", "CGT" = "B", "AGT" = "D",
                        "ACT" = "H", "ACG" = "V", "ACGT" = "N",
                        "N"
                    )
                }
                sequence <- gsub(token, iupac, sequence, fixed = TRUE)
            }
        }
        toupper(sequence)
    }

    write_fasta_format <- function(mat, filename, data_format = "dna") {
        sequences <- apply(mat, 1, paste, collapse = "")
        sequences <- vapply(
            sequences,
            function(sequence) fix_ambiguity_for_fasta(sequence, data_format),
            character(1)
        )
        out_lines <- as.vector(rbind(paste0(">", rownames(mat)), sequences))
        writeLines(out_lines, filename)
    }

    read_tnt_manual <- function(tnt_file) {
        lines <- readLines(tnt_file)
        xread_idx <- which(grepl("^xread", lines, ignore.case = TRUE))[1]
        if (is.na(xread_idx)) stop("Could not find 'xread' block in TNT file: ", tnt_file)
        lines <- lines[xread_idx:length(lines)]

        header_parts <- strsplit(trimws(lines[2]), "\\s+")[[1]]
        nchars <- as.numeric(header_parts[1])
        ntaxa <- as.numeric(header_parts[2])

        data_lines <- trimws(lines[-c(1, 2)])
        data_lines <- data_lines[nzchar(data_lines)]
        data_lines <- data_lines[!grepl("^&\\s*\\[", data_lines)]
        end_idx <- which(data_lines == ";")[1]
        if (!is.na(end_idx)) data_lines <- data_lines[seq_len(end_idx - 1)]

        seq_lines <- data_lines[seq_len(ntaxa)]
        seq_parts <- regmatches(seq_lines, regexec("^([^\\s]+)\\s+(.+)$", seq_lines, perl = TRUE))

        taxa_names <- vapply(seq_parts, `[`, character(1), 2)
        sequences <- vapply(seq_parts, `[`, character(1), 3)
        sequences <- gsub("\\s+", "", sequences, perl = TRUE)

        if (!any(grepl("[\\[\\{\\(]", sequences, perl = TRUE))) {
            parsed_nchars <- max(nchar(sequences, type = "chars"))
            mat <- matrix(unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE),
                nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE
            )
            rownames(mat) <- taxa_names
            return(mat)
        }

        token_pattern <- "\\[[^]]+\\]|\\{[^}]+\\}|\\([^)]*\\)|."
        mat_list <- regmatches(sequences, gregexpr(token_pattern, sequences, perl = TRUE))
        parsed_nchars <- max(lengths(mat_list))
        mat <- matrix(unlist(mat_list, use.names = FALSE), nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
        rownames(mat) <- taxa_names
        return(mat)
    }

    validate_output_files <- function(filenames) {
        output_dirs <- unique(dirname(filenames))
        output_dirs <- output_dirs[output_dirs != "."]
        missing_dirs <- output_dirs[!dir.exists(output_dirs)]
        if (length(missing_dirs) > 0) {
            stop("Output directory does not exist: ", paste(missing_dirs, collapse = ", "), call. = FALSE)
        }
    }

    tnt_fix_ambiguity <- function(sequence) {
        gsub("\\(([^)]+)\\)", "[\\1]", sequence, perl = TRUE)
    }

    write_tnt_lines <- function(lines, filename) {
        if (length(lines) > 0 && lines[1] == "nstates 32") {
            text <- paste0(lines[1], "\n", paste(lines[-1], collapse = "\r\n"), "\r\n")
            writeBin(charToRaw(text), filename)
        } else {
            writeLines(lines, filename, sep = "\r\n")
        }
    }

    write_tnt_format <- function(mat, filename, format = "standard") {
        ntaxa <- nrow(mat)
        nchars <- ncol(mat)
        lines <- c("xread", paste(nchars, ntaxa))
        if (format == "protein") lines <- c("nstates 32", lines)
        if (format == "protein") {
            lines <- c(lines, "", "& [prot]")
        } else if (format == "dna") {
            lines <- c(lines, "", "& [dna]")
        } else {
            lines <- c(lines, "", "& [num]")
        }

        for (i in 1:ntaxa) {
            taxon <- rownames(mat)[i]
            sequence <- tnt_fix_ambiguity(paste(mat[i, ], collapse = ""))
            lines <- c(lines, paste0(taxon, "\t", sequence))
        }
        lines <- c(lines, ";", "", "", "proc /;", "comments 0", ";", "", "")
        write_tnt_lines(lines, filename)
    }

    # --- Parse the list of taxa ---
    if (is.character(list) && length(list) == 1 && file.exists(list)) {
        taxa_list <- readLines(list)
    } else {
        taxa_list <- as.character(unlist(list))
    }
    taxa_list <- trimws(taxa_list)
    taxa_list <- taxa_list[nzchar(taxa_list)]

    if (length(taxa_list) == 0) stop("The provided taxa list is empty.")

    # --- Load input matrix ---
    data_format <- "unknown"
    if (is.matrix(input) || is.data.frame(input)) {
        mat <- as.matrix(input)
    } else {
        if (is.null(input_format)) {
            if (is.character(input)) {
                if (grepl("\\.tnt$", input, ignore.case = TRUE)) {
                    input_format <- "tnt"
                } else if (grepl("\\.nex(us)?$", input, ignore.case = TRUE)) {
                    input_format <- "nexus"
                } else if (grepl("\\.fa(s|sta)?$", input, ignore.case = TRUE)) {
                    input_format <- "fasta"
                } else {
                    stop("Could not auto-detect input_format. Provide input_format = 'nexus', 'tnt', or 'fasta'")
                }
            } else {
                stop("input must be a matrix, data frame, or file path.")
            }
        }

        if (input_format == "nexus") {
            data_format <- detect_format_from_nexus(input)
            mat <- TreeTools::ReadCharacters(input)
        } else if (input_format == "tnt") {
            data_format <- detect_format_from_tnt(input)
            mat <- read_tnt_manual(input)
        } else if (input_format == "fasta") {
            data_format <- detect_format_from_fasta(input)
            mat <- fasta_to_matrix(input, data_format)
        } else {
            stop("input_format must be 'nexus', 'tnt', or 'fasta'")
        }
        mat <- as.matrix(mat)
    }

    if (data_format == "unknown") {
        data_format <- detect_format_from_matrix(mat)
    }

    mat_taxa <- rownames(mat)
    if (is.null(mat_taxa)) stop("Input matrix must have row names (taxa names).")

    # --- Filter Logic ---
    if (filter == "keep") {
        selected_taxa <- intersect(mat_taxa, taxa_list)
    } else if (filter == "delete") {
        selected_taxa <- setdiff(mat_taxa, taxa_list)
    }

    if (length(selected_taxa) == 0) stop("Filtering removed all taxa from the matrix.")

    mat_filtered <- mat[selected_taxa, , drop = FALSE]
    removed_taxa <- setdiff(mat_taxa, selected_taxa)

    cat("Total taxa in original matrix:", length(mat_taxa), "\n")
    cat("Taxa in provided list:", length(taxa_list), "\n")
    cat("Filter mode:", filter, "\n")
    cat("Taxa removed after filter:", length(removed_taxa), "\n")
    if (length(removed_taxa) > 0) {
        cat("Removed taxa:\n")
        cat(paste0("  - ", removed_taxa, collapse = "\n"), "\n", sep = "")
    }
    cat("Taxa remaining after filter:", length(selected_taxa), "\n\n")

    # --- Output Logic ---
    if (!is.null(output)) {
        if (is.null(output_format)) {
            output_file_format <- if (grepl("\\.tnt$", output, ignore.case = TRUE)) {
                "tnt"
            } else if (grepl("\\.fa(s|sta)?$", output, ignore.case = TRUE)) {
                "fasta"
            } else {
                "nexus"
            }
        } else {
            output_file_format <- output_format
        }
        is_tnt_output <- output_file_format == "tnt"
        is_fasta_output <- output_file_format == "fasta"

        if (is_tnt_output && !grepl("\\.tnt$", output, ignore.case = TRUE)) {
            output <- paste0(output, ".tnt")
        } else if (is_fasta_output && !grepl("\\.fa(s|sta)?$", output, ignore.case = TRUE)) {
            output <- paste0(output, ".fasta")
        } else if (!is_tnt_output && !is_fasta_output && !grepl("\\.nex(us)?$", output, ignore.case = TRUE)) {
            output <- paste0(output, ".nexus")
        }

        validate_output_files(output)

        if (is_tnt_output) {
            write_tnt_format(mat_filtered, output, format = data_format)
            cat("Filtered matrix written to:", output, "\n")
            cat("Format: TNT (", data_format, ")\n")
        } else if (is_fasta_output) {
            if (!(data_format %in% c("dna", "protein"))) {
                stop("FASTA output is only supported for DNA or protein data.")
            }
            write_fasta_format(mat_filtered, output, data_format = data_format)
            cat("Filtered matrix written to:", output, "\n")
            cat("Format: FASTA (", data_format, ")\n")
        } else {
            nexus_data_format <- if (data_format %in% c("protein", "dna")) data_format else "standard"
            ape::write.nexus.data(mat_filtered, file = output, format = nexus_data_format, interleaved = FALSE)

            # Clean up the output like in filterSharedTaxa
            temp <- gsub("INTERLEAVE=NO", "", readLines(output))
            temp <- gsub("write.nexus.data.R", "RNODE", temp)
            writeLines(temp, output)

            cat("Filtered matrix written to:", output, "\n")
            cat("Format: NEXUS (", nexus_data_format, ")\n")
        }
        return(invisible(mat_filtered))
    } else {
        return(mat_filtered)
    }
}
