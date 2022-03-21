#!/usr/bin/Rscript --vanilla

args = commandArgs(trailingOnly=TRUE) #take in arguments for the script

het_file <- paste(args[1], ".het", sep="") #create variable for the het file
het_table <- read.table(het_file, header=T) # Read in the .het file

for (i in 1:nrow(het_table)) {
	het_table$temp[i] <- het_table$N.NM.[i] - het_table$O.HOM.[i]
}

m <- mean(het_table$temp) # Calculate the mean  
s <- sd(het_table$temp) # Calculate the SD

valid_by_number <- subset(het_table, temp <= m+3*s & temp >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean

for (i in 1:nrow(het_table)) {
	het_table$temp_2[i] <- het_table$temp[i] / het_table$N.NM.[i]
}

m2 <- mean(het_table$temp_2) # Calculate the second mean  
s2 <- sd(het_table$temp_2) # Calculate the second SD

valid_by_percent <- subset(het_table, temp_2 <= m2+3*s2 & temp_2 >= m2-3*s2) # Get any samples with F coefficient within 3 SD of the second mean

valid = merge(x=valid_by_number,y=valid_by_percent,by="IID") #make a new df with only samples that were valid by number and by percent

output_name <- paste(args[1], ".valid.sample", sep="") #add extension to output file variable

write.table(valid[,c(1,2)], output_name, quote=F, row.names=F) # print FID and IID for valid samples
