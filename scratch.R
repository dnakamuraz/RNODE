source("R/nexus2tnt.R")
source("R/tnt2nexus.R")

# Test nexus2tnt
nexus2tnt("testdata/048_MORPH_data.nex", "testdata/048_converted.tnt")
lines_tnt <- readLines("testdata/048_converted.tnt")
print("TNT snippet:")
print(head(lines_tnt, 10))
print(tail(lines_tnt, 5))

# Test tnt2nexus
tnt2nexus("testdata/059_MORPH_data.tnt", "testdata/059_converted.nexus")
lines_nex <- readLines("testdata/059_converted.nexus")
print("NEXUS snippet:")
print(head(lines_nex, 10))
print(tail(lines_nex, 5))
