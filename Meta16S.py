#!/usr/bin/python
# Copyright (C) 2018 Ion Torrent Systems, Inc. All Rights Reserved
#
# Plugin: Meta16S
# This plugin is developed for 16S Metagenomics data analysis
#
# Author: Lucius Zheng
# Last modified: 2018/11/07
#
# Main revisions:
# 	1. Add "Barcode Sample Settings +" in the plugin configure page to select which barcodes to process;
#	2. Update "Advanced options", as well as the default values for some parameters;
#	3. Include an HTML file named plan.html to enable plan configuration;
#	4. Add one download option for downloading all result files;
#	5. Add exceptions process to ensure a successful run.
#
# Last modified: 2019/02/14


import json
import os
from django.utils.functional import cached_property
from ion.plugin import *
import subprocess
from subprocess import check_output

from django.conf import settings
from django.template.loader import render_to_string

def createReport(reportName,reportTemplate,reportData):
	with open(reportName,'w') as bcsum:
		bcsum.write( render_to_string(reportTemplate,reportData) )

class Meta16S(IonPlugin):
	# The version number for this plugin
	version = "2.0.0.1"

	# this plugin can run on fullchip runs, thumbnail runs, and composite (merged via project page) runs
	# note that when the plugin is manually launched, only the 'launch' method will be called
	runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]

	# specify when the plugin is called.  For log parsing, stay simple and just get called when the run completes.
	# but can also be called before the run starts, at the block level, or after all other default plugins run
	runlevels = [RunLevel.DEFAULT]
	
	# a simple cached version of the start plugin property
	@cached_property
	def startplugin_json(self):
		return self.startplugin

	@cached_property
	def barcodes_json(self):
		with open('barcodes.json', 'r') as barcodes_handle:
			return json.load(barcodes_handle)
	
	def launch(self, data=None):
		"""This is the primary launch method for the plugin."""
	
		# configure django to use the templates folder		
		#settings.configure(TEMPLATE_DIRS=(self.startplugin["runinfo"]["plugin_dir"] + '/templates'),)
		
		if not settings.configured:
			settings.configure( DEBUG=False, TEMPLATE_DEBUG=False,
								INSTALLED_APPS=('django.contrib.humanize',),
								TEMPLATE_DIRS=(os.path.join(self.startplugin["runinfo"]["plugin_dir"],'templates'),) 
							)
		
		# define mothur command related environment variables
		mothur = self.startplugin_json['runinfo']['plugin_dir'] + "/mothur.1.39.5/mothur/mothur"
		mothur_bin = os.path.join(self.startplugin_json['runinfo']['plugin_dir'],'mothur.1.39.5/mothur')
		#mothur_bin = "/results/plugins/Meta16S_V3/scripts/mothur.1.39.5/mothur"
		database_dir = os.path.join(self.startplugin_json['runinfo']['plugin_dir'],'database')
		#database_dir = "/results/plugins/Meta16S_V3/scripts/database"
		script_dir = os.path.join(self.startplugin_json['runinfo']['plugin_dir'],'scripts')
		
		# save input parameters
		parameters_used = {}
		parameters_used["Number of reads"] = self.startplugin_json['pluginconfig']['num_of_reads']
		parameters_used["Forward primer"] = self.startplugin_json['pluginconfig']['primer_f']
		parameters_used["Reverse primer"] = self.startplugin_json['pluginconfig']['primer_r']
		parameters_used["Reference database"] = self.startplugin_json['pluginconfig']['database']
		parameters_used["Minimum average quality score"] = self.startplugin_json['pluginconfig']['qaverage']
		parameters_used["Maximum homopolymer length"] = self.startplugin_json['pluginconfig']['maxhomop']
		parameters_used["Maximum ambiguous bases"] = self.startplugin_json['pluginconfig']['maxambig']
		parameters_used["Minimum read length"] = self.startplugin_json['pluginconfig']['minlength']
		parameters_used["Maximum differences to primer"] = self.startplugin_json['pluginconfig']['pdiffs']
		parameters_used["Minimum reads required"] = self.startplugin_json['pluginconfig']['minnum']
		
		# save primer sequences info.
		primer_file= self.startplugin_json['runinfo']['results_dir'] +'/primers.txt'
		with open(primer_file,'w') as primer_handle:
			primer_f = "forward" + " " + self.startplugin_json['pluginconfig']['primer_f'] + "\n"
			primer_handle.write(primer_f)
			
			if not self.startplugin_json['pluginconfig']['primer_r'] == "":
				primer_r = "reverse" + " " + self.startplugin_json['pluginconfig']['primer_r'] + "\n"
				primer_handle.write(primer_r)
				
		# start to analyze bam files
		
		sample_num = 0
		for barcode_name, barcode_values in self.barcodes_json.iteritems():
			# do you work per barcode here!	
			
			# first check to see if the barcode was excluded using the frame work barcodes configuration table		
			selected = True
			barcodeData = self.startplugin_json['pluginconfig'].get('barcodetable',None)
			if barcodeData:
				#print(barcodeData)
				for bc in barcodeData:
					if  bc.get('barcode_name',"") == barcode_name:
						selected = bc.get('selected',True)
						break
			
			if not selected:
				continue
			
			print("Barcode Name: " + barcode_name)
			print("Bam Filepath: " + barcode_values['bam_filepath'])
			print("Read count: " + str(barcode_values['read_count']))
			
			if parameters_used.has_key("Barcodes selected"):
				parameters_used["Barcodes selected"].append(barcode_name)
			else:
				parameters_used["Barcodes selected"] = [barcode_name]
			
			# if no BAM file or file size is 0, then skip the sample
			if not os.path.exists(barcode_values['bam_filepath']):
				print "BAM file does not exist. We will skip the sample in the followed analysis.\n"
				continue
			if os.path.getsize(barcode_values['bam_filepath']) == 0:
				print "BAM file size is zero. We will skip the sample in the followed analysis.\n"
				continue
				
			# subsample reads to a fixed number, e.g., 10000
			
			if barcode_values['read_count'] < int(self.startplugin_json['pluginconfig']['minnum']):
				print "BAM file reads number is less than minimum required. We will skip the sample in the followed analysis.\n"
				continue
			
			# self.startplugin_json['pluginconfig']['num_of_reads'] type is unicode,
			# change to int type for comparison
			
			subsample_num = self.startplugin_json['pluginconfig']['num_of_reads']
			if barcode_values['read_count'] < int(self.startplugin_json['pluginconfig']['num_of_reads']):
				print "The selected number of analysis reads is larger than the raw reads count. We will use all reads in the followed analysis.\n"
				subsample_num = str(barcode_values['read_count'])
			
			arg1 = barcode_values['bam_filepath']
			cmd = "samtools view " + arg1 + " | shuf -n "+ subsample_num + "> output.sam"
			sampling_results = check_output(cmd, shell=True, cwd=self.startplugin_json['runinfo']['results_dir'])

			# convert sam/bam to fastq
			arg2 = barcode_name + '.fastq'
			fastq_results = check_output(["samtools", 'bam2fq', '-s', arg2, "output.sam"],cwd=self.startplugin_json['runinfo']['results_dir'])
			
			# read fastq file and create a fasta and quality file
			arg3 = "#fastq.info(fastq=" + arg2 + ")"
			mothur_read = check_output([mothur,arg3],cwd=self.startplugin_json['runinfo']['results_dir'])

			# remove the user-provided primers, barcodes, and sequences that drop below a quality threshold
			#
			# 	oligos: takes a file that can contain the sequences of the forward and reverse primers
			#	qaverage: tells mothur to calculate the average quality score for each sequence and 
			#		to remove those sequences that have an average below the value provided to the option
			#	flip=T to get the reverse complement of the sequences
			#	maxambig: cull those sequences that have ambiguous bases
			#	maxhomop: cap the homopolymer length
			#	minlength & maxlength: trim the sequence according their length
			# 	pdiffs: maximum number of differences to the primer sequence, default=0
			fasta_file = arg2.replace('.fastq','.fasta')
			qual_file = arg2.replace('.fastq','.qual')
			qaverage = self.startplugin_json['pluginconfig']['qaverage']
			maxhomop = self.startplugin_json['pluginconfig']['maxhomop']
			maxambig = self.startplugin_json['pluginconfig']['maxambig']
			minlength = self.startplugin_json['pluginconfig']['minlength']
			pdiffs = self.startplugin_json['pluginconfig']['pdiffs']
			
			arg4 = "#trim.seqs(fasta=" + fasta_file + ",oligos=" + primer_file + ",qfile=" + qual_file + ",flip=T,qaverage=" + qaverage + ",maxhomop=" + maxhomop + ",maxambig=" + maxambig + ",minlength=" + minlength + ",pdiffs=" + pdiffs + ",processors=4)"
			mothur_trim = check_output([mothur,arg4],cwd=self.startplugin_json['runinfo']['results_dir'])
						
			trim_fasta = fasta_file.replace('.fasta','.trim.fasta')
			
			# if no sequence, then skip this sample
			if os.path.getsize(self.startplugin_json['runinfo']['results_dir'] + "/" + trim_fasta) == 0:
				print "No sequence after trimming. We will skip the sample in the followed analysis.\n"
				continue
			
			sample_num += 1
			
			arg5 = "#make.group(fasta=" + trim_fasta + ", groups=" + barcode_name + ")"
			mothur_group = check_output([mothur,arg5],cwd=self.startplugin_json['runinfo']['results_dir'])
		
		# merge all groups and fasta
		merge_fasta = check_output(['cat *trim.fasta > all.fa'],shell=True,cwd=self.startplugin_json['runinfo']['results_dir'])
		merge_group = check_output(['cat *trim.groups > all.group'],shell=True,cwd=self.startplugin_json['runinfo']['results_dir'])
		
		# when all required input files are prepared, we run mothur pipeline to analyze the 16S rRNA data
		mothur_pipeline = os.path.join(self.startplugin_json['runinfo']['plugin_dir'],'mothur.sh')
		database = self.startplugin_json['pluginconfig']['database']
		
		print([mothur_pipeline,mothur_bin,database_dir,database,script_dir])
		
		try:
			mothur_16S = check_output([mothur_pipeline,mothur_bin,database_dir,database],cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		
		############################################################################################################################
		## Visualization analysis
		############################################################################################################################ 
		
		# Rarefaction plot for number of observed OTUs
		try:
			R_result = check_output(["Rscript", script_dir + "/rarePlot.R", 
									"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.rarefaction",
									"RarefactionCurve.png",
									"Number of Different OTUs"
									],
									cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		
		# Rarefaction plot for Shannon measure
		try:
			R_result = check_output(["Rscript", script_dir + "/rarePlot.R", 
									"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.r_shannon",
									"RarefactionCurve_shannon.png",
									"Rarefaction Measure: shannon"
									],
									cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		
		# Rarefaction plot for Chao1 measure
		try:
			R_result = check_output(["Rscript", script_dir + "/rarePlot.R", 
									"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.r_chao",
									"RarefactionCurve_chao.png",
									"Rarefaction Measure: chao1"
									],
									cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		
		# taxonomy binning		
		#
		# first, we need to reformat the outputs from mothur
		try:
			check_output(script_dir + "/prepare_input.sh",
						cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		#
		# then, get the relative abundance of each otu in a sample: Abundance / Total number of sequences in the group
		try:
			check_output(["Rscript", script_dir + "/relative_abundance.R"],
						cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		#
		# plot the distribution of taxonomic relative abundances across all taxonomic groups for all samples
		try:
			check_output(["Rscript", script_dir + "/taxonomic_binning.R"],
						cwd=self.startplugin_json['runinfo']['results_dir'])
		except subprocess.CalledProcessError as e:
			print e.output
		
		# beta diversity
		#
		# 1st: PCoA plot in 2-D
		
		if os.path.exists(self.startplugin_json['runinfo']['results_dir'] + "/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.pcoa.loadings"):
			try:
				check_output(["Rscript",script_dir + "/plotPCOA.R", 
						"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.pcoa.axes",
						"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.pcoa.loadings"],
						cwd=self.startplugin_json['runinfo']['results_dir'])
			except subprocess.CalledProcessError as e:
				print e.output

		
		# 2nd: PCoA plot in 3-D
		
		if os.path.exists(self.startplugin_json['runinfo']['results_dir'] + "/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.pcoa.loadings"):	
			try:
				check_output(["Rscript",script_dir + "/plotPCOA3d.R",
						"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.pcoa.axes",
						"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.pcoa.loadings"],
						cwd=self.startplugin_json['runinfo']['results_dir'])
			except subprocess.CalledProcessError as e:
				print e.output
		
		# Krona visualization
		if os.path.exists(self.startplugin_json['runinfo']['results_dir'] + "/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary"):		
			try:
				check_output(["python", script_dir + "/mothur_krona_XML.py",  
						"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary","output.xml"],
						cwd=self.startplugin_json['runinfo']['results_dir'])
			except subprocess.CalledProcessError as e:
				print e.output
		
			check_output("mkdir krona_analysis", shell=True, cwd=self.startplugin_json['runinfo']['results_dir'])
			check_output([script_dir + "/KronaTools-2.7/bin/ktImportXML","-o","krona_analysis/krona.html", "output.xml"], 
						cwd=self.startplugin_json['runinfo']['results_dir'])
					
		
		# visualizaiton of tree of samples
		if os.path.exists(self.startplugin_json['runinfo']['results_dir'] + "/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.tre"):		
			check_output(["./bin/figtree","-graphic","PNG",
					self.startplugin_json['runinfo']['results_dir'] + "/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.tre",
					self.startplugin_json['runinfo']['results_dir'] + "/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.tre.png"], 
					cwd= script_dir + "/FigTree_v1.4.3")
		
		###################################################################################################################
		## output in HTML
		###################################################################################################################
		
		adiversityData = []
		with open(self.startplugin['runinfo']['results_dir'] + '/all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.summary', 'r') as file_handle:
		#with open('all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.summary', 'r') as file_handle:
			next(file_handle)
			for line in file_handle:
				data_entry = {}
				data = line.strip().split("\t")
				data_entry['label']          = data[0]
				data_entry['group']          = data[1]
				data_entry['coverage']       = data[2]
				data_entry['nseqs']          = data[3]
				data_entry['sobs']           = data[4]
				data_entry['chao']           = data[5]
				data_entry['chao_lci']       = data[6]
				data_entry['chao_hci']       = data[7]
				data_entry['shannon']        = data[8]
				data_entry['shannon_lci']    = data[9]
				data_entry['shannon_hci']    = data[10]
				data_entry['invsimpson']     = data[11]
				data_entry['invsimpson_lci'] = data[12]
				data_entry['invsimpson_hci'] = data[13]
				
				adiversityData.append(data_entry)
		
		if database == "2015RDP": genome_name = "trainset14_032015.rdp.fasta"
		else: genome_name = "gg_13_8_99.fasta"
		render_context = {
			"autorefresh" : False,
			"genome_name" : genome_name,
			"library_type" : "16S fusion primers amplification",
			"reads_num": self.startplugin_json['pluginconfig']['num_of_reads'],
			"adiversityData" : adiversityData
		}
		
		cp_cmd = "cp " + os.path.join(self.startplugin_json['runinfo']['plugin_dir'],'templates/Meta16S_workflow.png') + " ./"
		check_output(cp_cmd,shell=True,cwd=self.startplugin_json['runinfo']['results_dir'])

		# save parameters into a json file
		json_file = self.startplugin_json['runinfo']['results_dir'] +'/parameters.json'
		with open(json_file,"w") as f:
			f.write(json.dumps(parameters_used))
		
		if sample_num > 1: 
			createReport(os.path.join(self.startplugin['runinfo']['results_dir'],'Meta16S_report.html'), 'barcode_summary_all.html', render_context )
			zip_cmd = "zip -r results.zip " + "*.png Meta16S_report.html parameters.json " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy " + \
					"all.OTU.summary.taxonomy " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.summary " + \
					"Taxonomic-Binning/ " + \
					"krona_analysis/ " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.tre"
			
		else:
			createReport(os.path.join(self.startplugin['runinfo']['results_dir'],'Meta16S_report.html'), 'barcode_specific.html', render_context )
			zip_cmd = "zip -r results.zip " + "*.png Meta16S_report.html parameters.json " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy " + \
					"all.OTU.summary.taxonomy " + \
					"all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.summary " + \
					"Taxonomic-Binning/ " + \
					"krona_analysis/"
		
		check_output(zip_cmd,shell=True,cwd=self.startplugin_json['runinfo']['results_dir'])
		results = '<h4> <a href="' + 'results.zip' + '" target="_blank">Download all result files</a></h4>'
		with open("download_block.html","w") as html_dl:
			html_dl.write("<html><body><pre>")
			html_dl.write(results)
			html_dl.write("</pre></body></html>")
		
		return True
	
	# Return list of columns you want the plugin table UI to show.
	# Columns will be displayed in the order listed.
	def barcodetable_columns(self):
		return [
			{ "field": "selected", "editable": True },
			{ "field": "barcode_name", "editable": False },
			{ "field": "sample", "editable": False } ]

# Devel use - running directly
if __name__ == "__main__":
	PluginCLI()
