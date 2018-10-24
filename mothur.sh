#!/bin/bash

BIN=$1
DATABASE_DIR=$2
DATABASE=$3

## basic sequence clean-up

# output: all.names,all.unique.fa
${BIN}/mothur "#unique.seqs(fasta=all.fa)"

# Output File Names: all.unique.align, all.unique.align.report, all.unique.flip.accnos
${BIN}/mothur "#align.seqs(fasta=all.unique.fa,reference=${DATABASE_DIR}/silva.bacteria.fasta,flip=T,gapopen=-1, processors=10)"

# Output File Names: all.unique.good.align, all.unique.bad.accnos, all.good.names, all.good.group
${BIN}/mothur "#screen.seqs(fasta=all.unique.align,optimize=start-end,criteria=90,group=all.group,name=all.names)"

# Output File Names: all.filter, all.unique.good.filter.fasta
${BIN}/mothur "#filter.seqs(fasta=all.unique.good.align,vertical=T,processors=10)"

# Output File Names: all.unique.good.filter.names,all.unique.good.filter.unique.fasta
${BIN}/mothur "#unique.seqs(fasta=all.unique.good.filter.fasta,name=all.good.names)"

# Output File Names:
# all.unique.good.filter.unique.precluster.fasta
# all.unique.good.filter.unique.precluster.names
# all.unique.good.filter.unique.precluster.IonXpress_004.map
# all.unique.good.filter.unique.precluster.IonXpress_005.map
# all.unique.good.filter.unique.precluster.IonXpress_006.map
# all.unique.good.filter.unique.precluster.IonXpress_007.map
# all.unique.good.filter.unique.precluster.IonXpress_008.map
# all.unique.good.filter.unique.precluster.IonXpress_009.map
# all.unique.good.filter.unique.precluster.IonXpress_010.map
# all.unique.good.filter.unique.precluster.IonXpress_011.map
# all.unique.good.filter.unique.precluster.IonXpress_012.map
# all.unique.good.filter.unique.precluster.IonXpress_013.map
# all.unique.good.filter.unique.precluster.IonXpress_014.map
# all.unique.good.filter.unique.precluster.IonXpress_015.map
${BIN}/mothur "#pre.cluster(fasta=all.unique.good.filter.unique.fasta,name=all.unique.good.filter.names,group=all.good.group,diffs=2,processors=10)"

# Output File Names:
# all.unique.good.filter.unique.precluster.denovo.uchime.chimeras
# all.unique.good.filter.unique.precluster.denovo.uchime.accn
${BIN}/mothur "#chimera.uchime(fasta=all.unique.good.filter.unique.precluster.fasta,
							name=all.unique.good.filter.unique.precluster.names,
							group=all.good.group,
							processors=10)"

# Output File Names:
# all.unique.good.filter.unique.precluster.pick.names
# all.unique.good.filter.unique.precluster.pick.fasta
# all.good.pick.group
${BIN}/mothur "#remove.seqs(accnos=all.unique.good.filter.unique.precluster.denovo.uchime.accnos,
							fasta=all.unique.good.filter.unique.precluster.fasta,
							name=all.unique.good.filter.unique.precluster.names,
							group=all.good.group)"

# classify your sequences using the Bayesian classifier based configured database
#
# We will also specific a cutoff of 80. This roughly means that we are 80%+ confident
# in our classification, a good thing if we’re going to remove groups. We want to be sure of what we’re removing.

if [ "${DATABASE}" = "2015RDP" ]
then
	# Output File Names:
	# all.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy
	# all.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
	# all.unique.good.filter.unique.precluster.pick.rdp.wang.flip.accnos
	${BIN}/mothur "#classify.seqs(fasta=all.unique.good.filter.unique.precluster.pick.fasta,
							name=all.unique.good.filter.unique.precluster.pick.names,
							group=all.good.pick.group,template=${DATABASE_DIR}/trainset14_032015.rdp.fasta,
							taxonomy=${DATABASE_DIR}/trainset14_032015.rdp.tax,
							cutoff=80,
							processors=10)"
	
	# Now that everything is classified, we want to remove our undesirables. 
    # We’ll just remove non-Bacteria domains and completely unclassifieds (“unknown”)
	#
	# Output File Names:
	# all.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy
	# all.unique.good.filter.unique.precluster.pick.pick.names
	# all.unique.good.filter.unique.precluster.pick.pick.fasta
	# all.good.pick.pick.group
	${BIN}/mothur "#remove.lineage(fasta=all.unique.good.filter.unique.precluster.pick.fasta,
							name=all.unique.good.filter.unique.precluster.pick.names,
							group=all.good.pick.group,
							taxonomy=all.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy,
							taxon=Mitochondria-Cyanobacteria_Chloroplast-unknown)"

	## define OTUs
	## Now we can define operational taxonomic units (OTUs), the sequencing proxy for a microbial species
						
	# firstly calculate distances between sequences (how difference they are from each other)
	# cutoff=0.15: tell mothur to not save any distances larger than 0.03
	# Output File Names: all.unique.good.filter.unique.precluster.pick.pick.dist
	#${BIN}/mothur "#dist.seqs(fasta=all.unique.good.filter.unique.precluster.pick.pick.fasta,cutoff=0.15,processors=8)"

	# and then clustering these distances based on a difference cutoff with cluster.split. In general, cutoffs are
	# 0.03 = 3% different or 97% similar ~ species
	# 0.05 = 5% different or 95% similar ~ genus
	# 0.1 = 10% different or 90% similar ~ family

	# cluster.split: use the taxonomic information to split the sequences into bins and then cluster within each bin. 
	# The Schloss lab have published results showing that if you split at the level of Order or Family, and cluster to a 0.03 cutoff, 
	# you’ll get just as good of clustering as you would with the "traditional" approach.
	${BIN}/mothur "#cluster.split(fasta=all.unique.good.filter.unique.precluster.pick.pick.fasta,
							name=all.unique.good.filter.unique.precluster.pick.pick.names,
							taxonomy=all.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,
							splitmethod=classify,
							taxlevel=4,
							cutoff=0.03,
							processors=10)"

	# generate OTU matrix
	# make.shared to combine these data with sample names to create a table we can understand. 
	# We will have mothur only give us the species-level OTUs label=0.03 
	# but you could ask for any level that you like (as long as it’s below the cutoff you may have used in dist.seqs)
	${BIN}/mothur "#make.shared(list=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,
							group=all.good.pick.pick.group,
							label=0.03)"


	# we probably also want to know the taxonomy for each of our OTUs: classify OTUs
	${BIN}/mothur "#classify.otu(list=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,
							group=all.good.pick.pick.group,
							name=all.unique.good.filter.unique.precluster.pick.names,
							taxonomy=all.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy,
							label=0.03)"
	
else
	${BIN}/mothur "#classify.seqs(fasta=all.unique.good.filter.unique.precluster.pick.fasta,
							name=all.unique.good.filter.unique.precluster.pick.names,
							group=all.good.pick.group,template=${DATABASE_DIR}/gg_13_8_99.fasta,
							taxonomy=${DATABASE_DIR}/gg_13_8_99.gg.tax,
							cutoff=80,
							processors=10)"
							
	${BIN}/mothur "#remove.lineage(fasta=all.unique.good.filter.unique.precluster.pick.fasta,
							name=all.unique.good.filter.unique.precluster.pick.names,
							group=all.good.pick.group,
							taxonomy=all.unique.good.filter.unique.precluster.pick.gg.wang.taxonomy,
							taxon=Mitochondria-Cyanobacteria_Chloroplast-unknown)"
	
	#${BIN}/mothur "#dist.seqs(fasta=all.unique.good.filter.unique.precluster.pick.pick.fasta,
	#					cutoff=0.15,processors=8)"

	${BIN}/mothur "#cluster.split(fasta=all.unique.good.filter.unique.precluster.pick.pick.fasta,
							name=all.unique.good.filter.unique.precluster.pick.pick.names,
							taxonomy=all.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy,
							splitmethod=classify,
							taxlevel=4,
							cutoff=0.03,
							processors=10)"

	${BIN}/mothur "#make.shared(list=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,
							group=all.good.pick.pick.group,
							label=0.03)"

	${BIN}/mothur "#classify.otu(list=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,
							group=all.good.pick.pick.group,
							name=all.unique.good.filter.unique.precluster.pick.names,
							taxonomy=all.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy,
							label=0.03)"
fi

## OTU based analysis

# A rarefaction curve plots the number of species as a function of the number of individuals sampled. 
# The curve usually begins with a steep slope, which at some point begins to flatten as fewer species 
# are being discovered per sample: the gentler the slope, the less contribution of the sampling to 
# the total number of operational taxonomic units or OTUs
#
# By default, the rarefaction.single() command uses 1,000 randomizations to generate the rarefaction 
# curve data for the observed number of OTUs (i.e. sobs).
#
# rarefaction curve data for observed number of OTUs (i.e. sobs) and other non-parametric richness estimators:
# shannon, chao, ace 
# output: all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.rarefaction,r_shannon,r_chao,r_ace
${BIN}/mothur "#rarefaction.single(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared,calc=sobs-shannon-chao-ace,processors=8)"


## alpha diversity: within sample diversity, to see overall diversity and richness of that community
#
# output: all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.summary
# richness (chao) & diversity (shannon)
${BIN}/mothur "#summary.single(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared,
						label=0.03,
						calc=coverage-nseqs-sobs-chao-shannon-invsimpson)"

## beta diversity: between sample diversity meaning every pairwise comparison of two samples has a unique value
#
# Generate a phylip-formatted distance matrix that describes the dissimilarity (1-similarity) among multiple group
# output: all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.dist
${BIN}/mothur "#dist.shared(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared,
							calc=braycurtis)"

# Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in 
# as few dimesnsions as possible
# output: two files ending in *pcoa.axes, and *pcoa.loadings. The *pcoa.loadings file will tell you 
#         what fraction of the total variance in the data are represented by each of the axes 
${BIN}/mothur "#pcoa(phylip=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.dist)"

# For principle component analysis (PCA), an overall microbial community is represented by a single point in an 
# xy- or xyz-plane. Two points that are closer to each other indicate that the overall microbiota of those two 
# samples are more similar than two points that are farther apart.
#${BIN}/mothur "# pca(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, label=0.03)"

# nMDS is similar to PCA in that it assigns a single xy(z) point to each sample relative to all other samples. 
# nMDS is more robust than PCA, however, because it permutes through calculations a number of times to find a 
# best fit. This fit is constrained within the number of axes you specify.

# Visualize beta-diversity metrics by heatmaps. Each box is a pairwise comparison between two
# samples. mothur has several beta-diversity calculator options. Here, we will use the Bray-Curtis metric.
# output: all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.braycurtis.heatmap.sim.svg
${BIN}/mothur "# heatmap.sim(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, label=0.03, calc=braycurtis)"

# Trees of samples
# the command will generate a newick-formatted tree file that describes the dissimilarity (1-similarity) among 
# multiple groups. Groups are clustered using the UPGMA algorithm using the distance between communities as 
# calculated using any of the calculators describing the similarity in community membership or structure
#
# output: all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.tre,
#         which can be visualized in software like TreeView or FigTree
#${BIN}/mothur "# tree.shared(phylip=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.dist)"
${BIN}/mothur "# tree.shared(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, label=0.03, calc=braycurtis)"

## OTU analyses
#
# Heatmaps of OTUs: each box represents the relative abundance of an OTU in a sample. You can scale the abundance
# to make differences and low abundance OTUs easier to see. We will use the default log10 scale and only visualize
# the 10 most abundant OTUs
# 
# output: all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.heatmap.bin.svg
${BIN}/mothur "# heatmap.bin(shared=all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, label=0.03, scale=log10, numotu=10)"

