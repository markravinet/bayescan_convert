# Genepop to BayesScan conversion
# Command line utility to use functions - must bundle with source file
source("./big_convert.R")
Args <- commandArgs(TRUE) 

infile <- Args[1]
outfile <- Args[2]

big_convert(infile = infile, outfile = outfile)