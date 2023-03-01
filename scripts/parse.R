
all_data <- readLines("./ex_energy.txt")
head(all_data)

# remove first two rows
all_data = all_data[-1]
all_data = all_data[-1]
head(all_data)

# remove long deliminator
all_data = gsub("\\*\\*","", all_data)
#head(all_data)

# remove first *
all_data = gsub("^\\*","", all_data)
#head(all_data)

# remove last *
all_data = gsub("\\*$","", all_data)
#head(all_data)

# change * comma *
all_data = gsub("\\*",",", all_data)
head(all_data)

# remove all empty spaces
all_data = all_data[which(all_data!="")]
head(all_data)

# remove all white space
all_data = gsub("\\s","",all_data)
head(all_data)

outfile <- file("ex_energy.csv")
writeLines(all_data,outfile)