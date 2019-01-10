REFS = data/references
FIGS = results/figures
TABLES = results/tables
PROC = data/process
FINAL = submission/

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'


################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analyses including the
# SILVA reference alignment and RDP reference taxonomy. Note that this code
# assumes that mothur is in your PATH. If not (e.g. it's in code/mothur/, you
# will need to replace `mothur` with `code/mothur/mothur` throughout the
# following code.
#
################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# 132 and described at http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/
# We will use the SEED v. 132, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N https://mothur.org/w/images/7/71/Silva.seed_v132.tgz
	tar xvzf Silva.seed_v132.tgz silva.seed_v132.align silva.seed_v132.tax
	mothur "#get.lineage(fasta=silva.seed_v132.align, taxonomy=silva.seed_v132.tax, taxon=Bacteria)"
	mv silva.seed_v132.pick.align $(REFS)/silva.seed.align
	mothur "#get.seqs(fasta=silva.seed_v132.align, accnos=$(REFS)/euks.accnos)"
	cat silva.seed_v132.pick.align >> $(REFS)/silva.seed.align
	rm Silva.seed_v132.tgz silva.seed_v132.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2017/03/15/RDP-v16-reference_files/

$(REFS)/trainset16_022016.% :
	wget -N https://www.mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz
	tar xvzf Trainset16_022016.pds.tgz trainset16_022016.pds
	mv trainset16_022016.pds/* $(REFS)/
	rm -rf trainset16_022016.pds
	rm Trainset16_022016.pds.tgz

# We need to get the Zymo mock community data; note that Zymo named the 5 operon of Salmonella twice
$(REFS)/zymo_mock.align : $(REFS)/silva.v4.align
	wget -N https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
	unzip ZymoBIOMICS.STD.refseq.v2.zip
	rm ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*itochondria_ssrRNA.fasta #V4 primers don't come close to annealing to these
	cat ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*fasta > zymo_temp.fasta
	sed '0,/Salmonella_enterica_16S_5/{s/Salmonella_enterica_16S_5/Salmonella_enterica_16S_7/}' zymo_temp.fasta > zymo.fasta
	mothur "#align.seqs(fasta=zymo.fasta, reference=data/references/silva.v4.align, processors=12)"
	mv zymo.align data/references/zymo_mock.align
	rm -rf zymo* ZymoBIOMICS.STD.refseq.v2* zymo_temp.fasta

################################################################################
#
# Part 2: Get and run data through mothur
#
#	Process fastq data through the generation of files that will be used in the
# overall analysis.
#
################################################################################

# here we take the raw fastq files for all of the files and process them through the generation of
# a shared file. we stop at different stages where files are needed for splitting off the mock
# community data from the stool data

# Run stool sequences from make.contigs through make.shared
data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared:\
			code/get_shared_stool.batch\
			data/references/silva.v4.align\
			data/references/trainset16_022016.pds.fasta\
			data/references/trainset16_022016.pds.tax\
			data/raw/stool.files
	mothur code/get_shared_stool.batch
	rm data/mothur/stool*map

data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.dist data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.ave-std.summary :\
		data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared\
		code/get_alpha_beta_stool.bash
	bash code/get_alpha_beta_stool.bash $^




# make.contigs; screen.seqs; unique; align (w/mock); filter; unique; classify.seqs; remove contaminants

data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.fasta data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table data/mothur/zymo_mock.filter.fasta:\
			code/get_good_mock.batch\
			data/references/silva.v4.align\
			$(REFS)/zymo_mock.align\
			data/references/trainset16_022016.pds.fasta\
			data/references/trainset16_022016.pds.tax\
			data/raw/mock.files
	mothur code/get_good_mock.batch


# error.seqs
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.error.summary:\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.fasta\
			data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table\
			data/mothur/zymo_mock.filter.fasta
	mothur "#seq.error(fasta=data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.fasta, count=data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table, reference=data/mothur/zymo_mock.filter.fasta, aligned=T)"


#need to remove those sequences that are more than 20 away from a reference
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.fasta data/mothur/mock.trim.contigs.good.unique.good.filter.pick.pick.count_table : \
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.fasta\
			data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.error.summary\
			code/remove_contaminants.R
	Rscript code/remove_contaminants.R


# precluster
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.fasta data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table:\
			code/run_precluster_mock.batch\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.fasta\
			data/mothur/mock.trim.contigs.good.unique.good.filter.pick.pick.count_table
	mothur code/run_precluster_mock.batch
	rm data/mothur/mock*map


# error.seqs
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.error.summary:\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.fasta\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table\
			data/mothur/zymo_mock.filter.fasta
	mothur "#seq.error(fasta=data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.fasta, count=data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table, reference=data/mothur/zymo_mock.filter.fasta, aligned=T)"


# chimera.vsearch
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.fasta data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.count_table:\
			code/run_vchime_mock.batch\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.fasta\
			data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table
	mothur code/run_vchime_mock.batch
	mv data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.pick.fasta data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.fasta
	mv data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.denovo.vsearch.pick.count_table data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.count_table


# dist.seqs; cluster based on chimera.vsearch output
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.shared :\
		code/get_vsearch_shared_mock.batch\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.fasta\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.count_table
	mothur code/get_vsearch_shared_mock.batch

data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.dist data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.merge.braycurtis.0.03.lt.ave.dist data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.merge.groups.ave-std.summary :\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.shared\
		code/get_alpha_beta_mock.bash
	bash code/get_alpha_beta_mock.bash $^


# perfect chimera removal
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.fasta data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.count_table : \
		code/perfect_chimera_removal.R\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.error.summary\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.fasta\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table
	Rscript code/perfect_chimera_removal.R


# dist.seqs; cluster based on perfect chimera removal; get shared file
data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.shared :\
		code/get_perfect_shared_mock.batch\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.fasta\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.count_table
	mothur code/get_perfect_shared_mock.batch

data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.braycurtis.0.03.lt.ave.dist data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.merge.braycurtis.0.03.lt.ave.dist data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.merge.groups.ave-std.summary :\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.shared\
		code/get_alpha_beta_mock.bash
	bash code/get_alpha_beta_mock.bash $^


# No sequencing errors, no error_chimera_rates
data/mothur/taxa_mapping.braycurtis.1.lt.ave.dist data/mothur/taxa_mapping.groups.ave-std.summary:\
		code/mock_perfect_sequencing.batch\
		data/mothur/taxa_mapping.shared
	bash code/mock_perfect_sequencing.batch


# No sequencing errors, no bias - through unique.seqs
data/mothur/zymo_mock.filter.pick.count_table data/mothur/zymo_mock.filter.pick.unique.fasta : \
		code/no_sequence_errors_unique.batch\
		data/mothur/zymo_mock.filter.fasta
	grep "18S" data/mothur/zymo_mock.filter.fasta | cut -c2- > data/mothur/18S.accnos
	mothur code/no_sequence_errors_unique.batch


data/mothur/zymo_mock.filter.pick.unique.precluster.opti_mcc.list data/mothur/zymo_mock.filter.pick.unique.precluster.count_table : \
		code/no_sequence_errors_cluster.batch\
		data/mothur/zymo_mock.filter.pick.count_table\
		data/mothur/zymo_mock.filter.pick.unique.fasta
	mothur code/no_sequence_errors_cluster.batch

data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes : code/get_pcoa.bash data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.dist
	bash $^

data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes : code/get_pcoa.bash data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.dist
	bash $^


################################################################################
#
# Part 3: Summary table generation
#
################################################################################

data/process/taxa_mapping.shared data/process/mock_bias.tsv data/process/mock_bias_salmonella.tsv:\
		code/bias_analysis.R\
		data/mothur/zymo_mock.filter.pick.unique.fasta\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.error.summary\
		data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table
	Rscript code/bias_analysis.R

data/process/error_chimera_rates.tsv : code/error_chimera_analysis.R\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.error.summary\
		data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.error.summary\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.count_table
	Rscript code/error_chimera_analysis.R

data/process/mock_alpha_diversity.tsv : code/mock_alpha.R\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.merge.groups.ave-std.summary\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.merge.groups.ave-std.summary\
		data/mothur/taxa_mapping.groups.ave-std.summary\
		data/mothur/zymo_mock.filter.pick.unique.precluster.opti_mcc.list\
		data/mothur/zymo_mock.filter.pick.unique.precluster.count_table
	Rscript code/mock_alpha.R

data/process/mock_beta_diversity.tsv : code/mock_beta.R\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.merge.braycurtis.0.03.lt.ave.dist\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.merge.braycurtis.0.03.lt.ave.dist
	Rscript code/mock_beta.R

data/process/mock_beta_drift.csv : code/mock_drift.R\
		data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.braycurtis.0.03.lt.ave.dist
	Rscript code/mock_drift.R

data/process/stool_alpha_diversity.csv : code/stool_alpha.R\
		data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.ave-std.summary
	Rscript code/stool_alpha.R

data/process/stool_chimera.csv : code/stool_chimera.R\
		data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.count_table\
		data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table
	Rscript code/stool_chimera.R

data/process/stool_beta_diversity.tsv : code/stool_beta.R\
		data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.dist
	Rscript code/stool_beta.R


################################################################################
#
# Part 6: Figure generation
#
################################################################################

results/figures/mock_error.pdf : code/plot_error_rate.R data/process/error_chimera_rates.tsv
	Rscript $<

results/figures/chimera_plots.pdf : code/plot_chimera_rate.R\
																		data/process/error_chimera_rates.tsv\
																		data/process/stool_chimera.tsv
	Rscript $<

results/figures/species_bias.pdf : code/plot_species_bias.R\
																		data/process/mock_bias.tsv
	Rscript $<

results/figures/salmonella_bias.pdf : code/plot_salmonella_bias.R\
																			data/process/mock_bias_salmonella.tsv
	Rscript $<

results/figures/mock_community.pdf : code/plot_mock_community.R\
				data/process/mock_alpha_diversity.tsv\
				data/process/mock_beta_diversity.tsv\
				data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes
	Rscript $<

results/figures/stool_community.pdf : code/plot_stool_community.R\
				data/process/stool_alpha_diversity.tsv\
				data/process/stool_beta_diversity.tsv\
				data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes
	Rscript $<

results/figures/drift.pdf : code/plot_drift.R\
				data/process/mock_beta_drift.csv\
				data/process/stool_beta_diversity.tsv
	Rscript $<


################################################################################
#
# Part 5: Pull it all together
#
# Render the manuscript
#
################################################################################

submission/figure_1.eps : results/figures/mock_error.pdf
	pdf2ps $< $@

submission/figure_2.eps : results/figures/chimera_plots.pdf
	pdf2ps $< $@

submission/figure_3.eps : results/figures/species_bias.pdf
	pdf2ps $< $@

submission/figure_4.eps : results/figures/mock_community.pdf
	pdf2ps $< $@

submission/figure_5.eps : results/figures/stool_community.pdf
	pdf2ps $< $@

submission/figure_6.eps : results/figures/drift.pdf
	pdf2ps $< $@


submission/figure_s1.eps : results/figures/salmonella_bias.pdf
	pdf2ps $< $@




$(FINAL)/manuscript.% : 			\ #include data files that are needed for paper don't leave this line with a : \
						$(FINAL)/mbio.csl\
						$(FINAL)/references.bib\
						$(FINAL)/manuscript.Rmd
	R -e 'render("$(FINAL)/manuscript.Rmd", clean=FALSE)'
	mv $(FINAL)/manuscript.knit.md submission/manuscript.md
	rm $(FINAL)/manuscript.utf8.md


write.paper : $(TABLES)/table_1.pdf $(TABLES)/table_2.pdf\ #customize to include
				$(FIGS)/figure_1.pdf $(FIGS)/figure_2.pdf\	# appropriate tables and
				$(FIGS)/figure_3.pdf $(FIGS)/figure_4.pdf\	# figures
				$(FINAL)/manuscript.Rmd $(FINAL)/manuscript.md\
				$(FINAL)/manuscript.tex $(FINAL)/manuscript.pdf
