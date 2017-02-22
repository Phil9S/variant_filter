rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

af <- as.character(args[1])

var <- as.character(args[2])
sample <- as.character(args[3])

print(af)
print(var)
print(sample)

aftable <- read.table(af, header = TRUE,stringsAsFactors = FALSE)



val <- aftable[var,sample]

print(valtrue)
print(val)
