#' @title renameTaxa
#' @name renameTaxa
#' @description Renames terminals (taxa) in a phylogenetic matrix (TNT, NEXUS, or an R matrix)
#'   using a provided synonym dictionary. It can optionally delete taxa that are not specified
#'   in the synonym file and output the resulting matrix to a local file.
#' @author [Your Name]
#' @param input A local file path (NEXUS or TNT) or a matrix/data.frame already loaded in R.
#' @param input_format Format of the input file: 'nexus' or 'tnt'. Detected automatically if NULL (default).
#' @param synonym_file Path to a text file containing the synonym list. In each row, the first string
#'   is the new name, and subsequent strings (separated by spaces or tabs) are the old names/synonyms.
#' @param delete_missing Logical. If TRUE, deletes all terminals in the input matrix that are not present
#'   in the `synonym_file`. Default is FALSE.
#' @param output Output file path (e.g. "output.tnt"). If NULL, the result is returned as an R object only.
#' @param output_format Output format: 'nexus' or 'tnt'. If NULL, inferred from the extension of `output`. Ignored when output is NULL.
#' @return Invisibly returns the matrix with renamed (and optionally filtered) terminals.
#' @examples
#' \dontrun{
#' # Rename taxa in a TNT file and write the result back to TNT
#' renameTaxa(
#'     input = "input.tnt",
#'     synonym_file = "synonyms.txt",
#'     delete_missing = TRUE,
#'     output = "renamed_output.tnt"
#' )
#'
#' # Read from a loaded R matrix and return an R object only
#' mat <- renameTaxa(
#'     input = my_r_matrix,
#'     synonym_file = "synonyms.txt"
#' )
#' }
#'
#' @export
renameTaxa <- function(input,
                       input_format = NULL,
                       synonym_file,
                       delete_missing = FALSE,
                       output = NULL,
                       output_format = NULL) {
    # ---------------------------------------------------------------------------
    # Internal helpers
    # ---------------------------------------------------------------------------

    # Auto-detect format from file extension
    detect_format_from_ext <- function(path) {
        if (grepl("\\.tnt$", path, ignore.case = TRUE)) "tnt" else "nexus"
    }

    # Detect data type (standard / dna / protein) from a NEXUS file
    detect_type_from_nexus <- function(nexus_file) {
        lines <- tolower(readLines(nexus_file))
        fmt <- grep("\\bformat\\b", lines, value = TRUE)
        if (length(fmt) > 0) {
            if (grepl("datatype\\s*=\\s*protein|\\bprotein\\b", fmt[1])) {
                return("protein")
            }
            if (grepl("datatype\\s*=\\s*dna|\\bdna\\b", fmt[1])) {
                return("dna")
            }
            if (grepl("standard", fmt[1])) {
                return("standard")
            }
        }
        "standard"
    }

    # Detect data type from a TNT file
    detect_type_from_tnt <- function(tnt_file) {
        lines <- readLines(tnt_file)
        lines_lc <- tolower(lines)
        if (any(grepl("nstates\\s+(prot|32)", lines_lc)) || any(grepl("&\\s*\\[prot\\]", lines_lc))) {
            return("protein")
        }
        if (any(grepl("nstates\\s+dna", lines_lc)) || any(grepl("&\\s*\\[dna\\]", lines_lc))) {
            return("dna")
        }
        "standard"
    }

    # Guess type if input is an R matrix (fallback)
    guess_type_from_mat <- function(m) {
        chars <- unique(as.character(m))
        chars <- chars[!chars %in% c("?", "-", "N", "n")]
        if (length(chars) == 0) {
            return("standard")
        }
        if (all(tolower(chars) %in% c("a", "c", "g", "t", "u", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v"))) {
            return("dna")
        }
        if (any(toupper(chars) %in% LETTERS)) {
            return("protein")
        }
        return("standard")
    }

    # Read a TNT file, returning a character matrix (taxa x characters)
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

        if (length(data_lines) < ntaxa) {
            stop("Expected ", ntaxa, " taxa but found ", length(data_lines), " sequence lines in: ", tnt_file)
        }
        seq_lines <- data_lines[seq_len(ntaxa)]

        seq_parts <- regmatches(seq_lines, regexec("^([^\\s]+)\\s+(.+)$", seq_lines, perl = TRUE))
        parsed <- vapply(seq_parts, length, integer(1)) == 3
        if (!all(parsed)) {
            stop("Could not parse taxon/sequence on line(s): ", paste(which(!parsed), collapse = ", "))
        }

        taxa_names <- vapply(seq_parts, `[`, character(1), 2)
        sequences <- vapply(seq_parts, `[`, character(1), 3)
        sequences <- gsub("\\s+", "", sequences, perl = TRUE)

        # Fast path: no polymorphic tokens
        if (!any(grepl("[\\[\\{\\(]", sequences, perl = TRUE))) {
            seq_widths <- nchar(sequences, type = "chars")
            parsed_nchars <- max(seq_widths)
            if (!all(seq_widths == parsed_nchars)) {
                split_sequences <- strsplit(sequences, "", fixed = TRUE)
                split_sequences <- lapply(split_sequences, function(x) {
                    if (length(x) < parsed_nchars) c(x, rep("?", parsed_nchars - length(x))) else x
                })
                mat <- matrix(unlist(split_sequences, use.names = FALSE), nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE)
                rownames(mat) <- taxa_names
                warning("TNT sequences have unequal lengths; shorter sequences padded with '?'.")
                return(mat)
            }
            mat <- matrix(unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE),
                nrow = ntaxa, ncol = parsed_nchars, byrow = TRUE
            )
            rownames(mat) <- taxa_names
            return(mat)
        }

        # Slow path: polymorphic / inapplicable tokens
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
        return(mat)
    }

    # Convert (01)-style ambiguity to TNT [01]
    tnt_fix_ambiguity <- function(sequence) {
        gsub("\\(([^)]+)\\)", "[\\1]", sequence, perl = TRUE)
    }

    # Write a single block TNT file formatted to mimic the user's template
    write_tnt_single <- function(mat, filename, data_type = "standard") {
        ntaxa <- nrow(mat)
        nchars <- ncol(mat)

        lines <- character()
        if (data_type == "dna") lines <- c("nstates dna; ")
        if (data_type == "protein") lines <- c("nstates prot; ")

        lines <- c(lines, paste("xread", nchars, ntaxa))

        for (i in seq_len(ntaxa)) {
            seq <- tnt_fix_ambiguity(paste(mat[i, ], collapse = ""))
            lines <- c(lines, paste0(rownames(mat)[i], "\t", seq, " "))
        }
        lines <- c(lines, "; ", "proc /; ", "comments 0 ;")
        writeLines(lines, filename)
    }

    # ---------------------------------------------------------------------------
    # Load input
    # ---------------------------------------------------------------------------
    data_type <- "standard"

    if (is.matrix(input) || is.data.frame(input)) {
        mat <- as.matrix(input)
        data_type <- guess_type_from_mat(mat)
    } else {
        if (!file.exists(input)) stop("Input file not found: ", input)

        if (is.null(input_format)) input_format <- detect_format_from_ext(input)
        input_format <- match.arg(tolower(input_format), c("nexus", "tnt"))

        if (input_format == "nexus") {
            data_type <- detect_type_from_nexus(input)
            mat <- as.matrix(TreeTools::ReadCharacters(input))
        } else {
            data_type <- detect_type_from_tnt(input)
            mat <- read_tnt_manual(input)
        }
    }

    if (is.null(rownames(mat))) {
        stop("The input matrix must have row names (taxon names).")
    }

    # ---------------------------------------------------------------------------
    # Load and parse Synonym Dictionary
    # ---------------------------------------------------------------------------
    if (!file.exists(synonym_file)) stop("Synonym file not found: ", synonym_file)
    syn_lines <- readLines(synonym_file)
    syn_lines <- syn_lines[nzchar(trimws(syn_lines))]

    # Dictionary mapping: old_name/synonym -> new_name
    taxa_dict <- character()
    for (line in syn_lines) {
        parts <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(parts) > 0) {
            new_name <- parts[1]
            # Map the new name to itself
            taxa_dict[new_name] <- new_name
            # Map all synonyms to the new name
            if (length(parts) > 1) {
                for (syn in parts[-1]) {
                    taxa_dict[syn] <- new_name
                }
            }
        }
    }

    # ---------------------------------------------------------------------------
    # Rename and Filter Terminals
    # ---------------------------------------------------------------------------
    old_names <- rownames(mat)
    new_names <- character(length(old_names))
    keep <- rep(TRUE, length(old_names))

    for (i in seq_along(old_names)) {
        nm <- old_names[i]
        if (nm %in% names(taxa_dict)) {
            new_names[i] <- taxa_dict[nm]
        } else {
            new_names[i] <- nm # keep old name if not found in dict
            if (delete_missing) {
                keep[i] <- FALSE
            }
        }
    }

    mat <- mat[keep, , drop = FALSE]
    rownames(mat) <- new_names[keep]

    if (nrow(mat) == 0) {
        warning("All taxa were deleted! Either the input matrix was empty or no taxa matched the synonym file (and delete_missing = TRUE).")
    } else {
        cat(sprintf("Original taxa: %d | Remaining taxa: %d\n\n", length(old_names), nrow(mat)))
    }

    # ---------------------------------------------------------------------------
    # Write Output File (if requested)
    # ---------------------------------------------------------------------------
    if (!is.null(output) && nrow(mat) > 0) {
        if (is.null(output_format)) {
            output_format <- if (grepl("\\.tnt$", output, ignore.case = TRUE)) "tnt" else "nexus"
        }
        output_format <- match.arg(tolower(output_format), c("nexus", "tnt"))

        out_dir <- dirname(output)
        if (out_dir != "." && !dir.exists(out_dir)) {
            stop("Output directory does not exist: ", out_dir)
        }

        if (output_format == "tnt") {
            if (!grepl("\\.tnt$", output, ignore.case = TRUE)) output <- paste0(output, ".tnt")
            write_tnt_single(mat, output, data_type)
            cat("Renamed matrix written to:", output, "\nFormat: TNT\n")
        } else {
            if (!grepl("\\.(nexus|nex)$", output, ignore.case = TRUE)) output <- paste0(output, ".nexus")
            ape::write.nexus.data(mat, file = output, format = data_type, interleaved = FALSE)

            # Clean up Ape's output styling to match cleaner logic
            tmp <- readLines(output)
            tmp <- gsub("INTERLEAVE=NO", "", tmp)
            tmp <- gsub("write.nexus.data.R", "RNODE", tmp)
            writeLines(tmp, output)
            cat("Renamed matrix written to:", output, "\nFormat: NEXUS (", data_type, ")\n")
        }
    }

    invisible(mat)
}
