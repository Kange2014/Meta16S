#' Script: Taxonomic Binning
#'
#' Provides an overview of sample-specific relative abundances for all taxonomic levels
#'
#' Input:
#'    Write the name of the examined OTU file
#'
#' Output: The script is generating seven tab-delimited files
#' 1. Relative taxonomic abundance at the kingdom level for each sample
#' 2. Relative taxonomic abundance at the phyla level for each sample
#' 3. Relative taxonomic abundance at the class level for each sample
#' 4. Relative taxonomic abundance at the order level for each sample
#' 5. Relative taxonomic abundance at the family level for each sample
#' 6. Relative taxonomic abundance at the genera level for each sample
#' 7. Relative taxonomic abundance at the species level for each sample
#' 8. Relative taxonomic abundance at all taxonomic levels for each sample
#' 
#' Graphical Output:
#' 9. Distribution of taxonomic relative abundances across all taxonomic groups for all samples
#'
#' Concept:
#' Taxonomic information is split into the different taxonomic levels for each OTU
#' Relative abundance values of OTUs belonging to the same taxonomic group are summed up
#' Individual taxonomic composition for each sample is generated
#'
#' Note:
#' If taxonomic information at a specific level is missing, the entry is replaced by
#' the last available taxonomic level including the prefix "unknown_"

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please give the file name of the OTU-table containing relative abundances and taxonomic classification 
otu_file <- "OTUs_Table-rel-tax.tab"

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ###### 
##################################################################################

###################            Read input table              ####################

# Load the tab-delimited file containing the abundances and taxonomic information to be checked (rownames in the first column)
otu_table <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Create a dataframe with a number of rows identical to the number of OTUs in the dataset
taxonomy <- otu_table[,dim(otu_table)[2]]

# Test if the taxonomy column is in the correct format (delimited by semicolon)
# IMPORTANT: taxomomy classification from different databases (2015RDP, Greengenes) shows different formats
if(any(grepl("(?:[^;]*;){6}", taxonomy))==FALSE) {

#Send error message if taxonomy is not in the right format
  stop("Wrong number of taxonomic classes\n

Taxonomic levels have to be separated by semicolons (six or seven in total). 
IMPORTANT: if taxonomic information at any level is missing, the semicolons are still needed:\n
       
      e.g.k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__guillouiae;
      e.g.k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;;
	  e.g.Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Klebsiella;")
} else { 

# Delete the taxonomy column from the OTU-table
#otuFile <- otu_table[,c(1:dim(otu_table)[2] - 1)]  ## this will generate a numeric vector when otu_table only has 2 columns (one analyzed sample)
otuFile <- otu_table
otuFile[,dim(otu_table)[2]] <- NULL  ## this ensures a data.frame variable for downstream analysis

# Initialize empty dataframe
taxonomy_new <- NULL

# Split taxonomic information in its taxonomic classes
# Kingdom - Phylum - Class - Family - Order - Genus - Species
taxonomy_new <- strsplit(x = as.character(taxonomy), ";")

# Adjust dataframe with modified taxonomic information
taxonomy_new <- t(as.data.frame(taxonomy_new))
row.names(taxonomy_new) <- row.names(otuFile)

# record the length of taxonomy classes 
taxonomy_len <- length(taxonomy_new[1,])

# Add level information to all taxonomies when reference database is RDP which contains 6 tax levels
if (taxonomy_len == 6){

	# For taxonomies related to kingdom level
	taxonomy_new[,1] <- sub("^","k__",taxonomy_new[,1])

	for (i in 1:length(taxonomy_new[,1])) {

		# Save the position where the first string matches "unclassified"
		value <- grep("unclassified",taxonomy_new[i,])[1]
  
		if(is.na(value)){
			# For taxonomies related to phylum level
			taxonomy_new[i,2] <- sub("^","p__",taxonomy_new[i,2])

			# For taxonomies related to class level
			taxonomy_new[i,3] <- sub("^","c__",taxonomy_new[i,3])

			# For taxonomies related to order level
			taxonomy_new[i,4] <- sub("^","o__",taxonomy_new[i,4])

			# For taxonomies related to family level
			taxonomy_new[i,5] <- sub("^","f__",taxonomy_new[i,5])

			# For taxonomies related to genus level
			taxonomy_new[i,6] <- sub("^","g__",taxonomy_new[i,6])
		}
  
		else{
			if (value == 2) {
				taxonomy_new[i,2:6] <- sub("^","k__",taxonomy_new[i,2:6])
			}
  
			if (value == 3) {
				# For taxonomies related to phylum level
				taxonomy_new[i,2:6] <- sub("^","p__",taxonomy_new[i,2:6])
			}
			
			if (value == 4) {
				# For taxonomies related to phylum level
				taxonomy_new[i,2] <- sub("^","p__",taxonomy_new[i,2])
	
				# For taxonomies related to class level
				taxonomy_new[i,3:6] <- sub("^","c__",taxonomy_new[i,3:6])
			}
  
			if (value == 5) {
				# For taxonomies related to phylum level
				taxonomy_new[i,2] <- sub("^","p__",taxonomy_new[i,2])
	
				# For taxonomies related to class level
				taxonomy_new[i,3] <- sub("^","c__",taxonomy_new[i,3])
	
				# For taxonomies related to order level
				taxonomy_new[i,4:6] <- sub("^","o__",taxonomy_new[i,4:6])
			}
  
			if (value == 6) {
				# For taxonomies related to phylum level
				taxonomy_new[i,2] <- sub("^","p__",taxonomy_new[i,2])
	
				# For taxonomies related to class level
				taxonomy_new[i,3] <- sub("^","c__",taxonomy_new[i,3])
	
				# For taxonomies related to order level
				taxonomy_new[i,4] <- sub("^","o__",taxonomy_new[i,4])
	
				# For taxonomies related to family level
				taxonomy_new[i,5:6] <- sub("^","f__",taxonomy_new[i,5:6])
			}
		}
	}
}

#################################################################################

# Create list with taxonomic information for each taxonomy level

if(taxonomy_len == 6){
	class_list <- list(
		unique(taxonomy_new[,1]),unique(taxonomy_new[,2]),
		unique(taxonomy_new[,3]),unique(taxonomy_new[,4]),
		unique(taxonomy_new[,5]),unique(taxonomy_new[,6])
	)
}else{
	class_list <- list(
		unique(taxonomy_new[,1]),unique(taxonomy_new[,2]),
		unique(taxonomy_new[,3]),unique(taxonomy_new[,4]),
		unique(taxonomy_new[,5]),unique(taxonomy_new[,6]),
		unique(taxonomy_new[,7])
	)
}

  
# Clone the created list for further processing
sample_list <- class_list
list_length <- NULL

# Iterate through all seven taxonomy levels
for (a in 1:taxonomy_len) {
  
  lis <- lapply(class_list[a], lapply, length)
  names(lis)<-lapply(class_list[a],length)
  
  # Individual number of taxonomies for each taxonomic level
  num_taxa <- as.integer(names(lis))
  list_length[a] <- num_taxa
  
  # Iterate through taxonomic class specific taxonomies
  for (b  in 1:num_taxa) {
    
    # Initialize list with the value zero for all taxonomies
    sample_list[[a]][[b]] <- list(rep.int(0,dim(otuFile)[2]))
    
  }
}

#################################################################################
#################################################################################
# Save relative abundances of all samples for each taxonomy

matrix <- data.matrix(otuFile)

# Iterate through all taxonomic levels
for (m in 1:taxonomy_len) {
	
	# Iterate through all taxonomic classes at a particular taxonomic level
	for(n in 1:length(class_list[[m]])){
		
		# All rows with a particular taxonomic class of n samples
		sub_sample_tax <-(subset(matrix,taxonomy_new[,m] == class_list[[m]][n]))
		
		# Calculate the summed up relative abundances for the particular taxonomic class for n samples
		sample_list[[m]][[n]] <- list(colSums(sub_sample_tax))
	}
}

#################################################################################
######                         Write output                                ######
#################################################################################

# Generate tables for each taxonomic class

##Kingdom table
# Create table with taxonomic information (kingdom level)
kingdom <-  matrix(unlist(sample_list[[1]]),nrow = dim(otuFile)[2],ncol = list_length[1],dimnames = list(names(otuFile),unlist(class_list[[1]])))
kingdom <- (t(kingdom))

##Phylum table
# Create table with taxonomic information (phylum level)
phyla <- matrix(unlist(sample_list[[2]]),nrow = dim(otuFile)[2],ncol = list_length[2],dimnames = list(names(otuFile),unlist(class_list[[2]])))
phyla <- (t(phyla))

# Order table according to taxonomic name (descending)
phyla <- phyla[order(row.names(phyla)),,drop=F] # drop=F to ensure a single column data.frame

## Class table
# Create table with taxonomic information (class level)
classes <- matrix(unlist(sample_list[[3]]), nrow = dim(otuFile)[2], ncol = list_length[3], dimnames = list(names(otuFile),unlist(class_list[[3]])))
classes <- (t(classes))

# Order dataframe according to taxonomic name (descending)
classes <- classes[order(row.names(classes)),,drop=F]

## Orders
# create table with taxonomic information (Order)
orders <-matrix(unlist(sample_list[[4]]),nrow = dim(otuFile)[2],ncol = list_length[4],dimnames = list(names(otuFile),unlist(class_list[[4]])))
orders <- (t(orders))

# Order dataframe according to taxonomic name (descending)
orders <- orders[order(row.names(orders)),,drop=F]

## Family table
# Create table with taxonomic information (family level)
families <-matrix(unlist(sample_list[[5]]),nrow = dim(otuFile)[2],ncol = list_length[5],dimnames = list(names(otuFile),unlist(class_list[[5]])))
families <- (t(families))

# Order dataframe according to taxonomic name (descending)
families <- families[order(row.names(families)),,drop=F]

## Genus level
# Create table with taxonomic information (generum level)
genera <- matrix(unlist(sample_list[[6]]),nrow = dim(otuFile)[2],ncol = list_length[6],dimnames = list(names(otuFile),unlist(class_list[[6]])))
genera <- (t(genera))

# Order dataframe according to taxonomic name (descending)
genera <- genera[order(row.names(genera)),,drop=F]


## Species level
# Create table with taxonomic information (species level)
if(taxonomy_len == 7) {
	species <- matrix(unlist(sample_list[[7]]),nrow = dim(otuFile)[2],ncol = list_length[7],dimnames = list(names(otuFile),unlist(class_list[[7]])))
	species <- (t(species))
	# Order dataframe according to taxonomic name (descending)
	species <- species[order(row.names(species)),,drop=F]
}

# Merge all dataframes
if(taxonomy_len == 7) {
	tax_summary <-rbind.data.frame(kingdom,phyla,classes,orders,families,genera,species)
}else{
	tax_summary <-rbind.data.frame(kingdom,phyla,classes,orders,families,genera)
}

# Identify duplicates and remove them
tax_summary <- tax_summary[!duplicated(row.names(tax_summary)),]

################################################################################
######                        Write Output Files                           ######
#################################################################################

# Create a directory 
dir.create("Taxonomic-Binning")

# Set path for all outputs to the new directory
setwd("Taxonomic-Binning")

# Write output files for taxonomic composition of every sample
write.table(kingdom,"0.Kingdom.all.tab",sep = "\t",col.names = NA)
write.table(phyla,"1.Phyla.all.tab",sep = "\t",col.names = NA)
write.table(classes,"2.Classes.all.tab",sep = "\t",col.names = NA)
write.table(orders,"3.Orders.all.tab",sep = "\t",col.names = NA)
write.table(families,"4.Families.all.tab",sep = "\t",col.names = NA)
write.table(genera,"5.Genera.all.tab",sep = "\t",col.names = NA)
if(taxonomy_len == 7) { write.table(species,"6.Species.all.tab",sep = "\t",col.names = NA) }
write.table(tax_summary,"tax.summary.all.tab",sep = "\t",col.names = NA)
suppressWarnings (try(write.table(tax_summary, "../../5.Serial-Group-Comparisons/tax.summary.all.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))

#################################################################################
######                        Write Graphical Output                       ######
#################################################################################

#Kingdom
png("taxonomic-kingdom.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))

#k_col=distinctColorPalette(dim(kingdom)[1])
k_col=rainbow(dim(kingdom)[1])
k_col=sample(k_col)
barplot(kingdom,col=k_col, cex.names=0.5, ylab="cumulative relative abundance (%)", las=2, main="Taxonomic binning at Kingdom level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(kingdom)),cex=0.7,col = rev(k_col),pch = 16,pt.cex = 1.2)

dev.off()

#Phyla
png("taxonomic-phyla.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))


#p_col=distinctColorPalette(dim(phyla)[1])
p_col=rainbow(dim(phyla)[1])
p_col=sample(p_col)
barplot(phyla,col=p_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Phyla level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(phyla)),cex=0.7,col = rev(p_col),pch = 16,pt.cex = 1.2)

dev.off()

#Classes
png("taxonomic-class.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))

c_col=rainbow(dim(classes)[1])
c_col=sample(c_col)
barplot(classes,col=c_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Class level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(classes)),cex=0.7,col = rev(c_col),pch = 16,pt.cex = 1.2)

dev.off()

#Orders
png("taxonomic-order.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))

o_col=rainbow(dim(orders)[1])
o_col=sample(o_col)
barplot(orders,col=o_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Order level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(orders)),cex=0.7,col = rev(o_col),pch = 16,pt.cex = 1.2)

dev.off()

#Families
png("taxonomic-family.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))

f_col=rainbow(dim(families)[1])
f_col=sample(f_col)
barplot(families,col=f_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Family level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(families)),cex=0.7,col = rev(f_col),pch = 16,pt.cex = 1.2)

dev.off()

#Genera
png("taxonomic-genus.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))

g_col=rainbow(dim(genera)[1])
g_col=sample(g_col)
barplot(genera,col=g_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Genus level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(genera)),cex=0.7,col = rev(g_col),pch = 16,pt.cex = 1.2)

dev.off()

if(taxonomy_len == 7) {

#Species
png("taxonomic-species.png",res=600,width=4800, height=4800)
par(xpd=T, mar=par()$mar+c(0,0,0,9))

g_col=rainbow(dim(species)[1])
g_col=sample(g_col)
barplot(species,col=g_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Species level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(species)),cex=0.7,col = rev(g_col),pch = 16,pt.cex = 1.2)

dev.off()
}

}

#################################################################################
######                           End of Script                             ######
#################################################################################
