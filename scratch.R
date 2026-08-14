# Test the parsing logic
input_constraint <- "(Pristimantis_fenestratus_MTR37091_MT,Ischnocnema_parva_MNRJ51356)"

if (is.character(input_constraint)) {
  input_constraint <- lapply(input_constraint, function(x) {
    x <- gsub("[()]", "", x)
    x <- strsplit(x, ",")[[1]]
    trimws(x)
  })
}
print(input_constraint)
