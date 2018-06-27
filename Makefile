M_SAMPLING = 50 100 500 1000 5000 10000
H_SAMPLING = 1000 5000 10000 15000 20000 
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
# SILVA reference alignment and RDP reference taxonomy.
#
################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# We will use the SEED v. 123, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz

################################################################################
#
# Part 2: Run data through mothur
#
#	Process fastq data through the generation of files that will be used in the
# overall analysis.
#
################################################################################

# Set up specific Taq variables
TAQ_USED = acc_data k_data phu_data pl_data q5_data
TAQ_USED_PATH=$(addprefix $(PROC)/,$(H_SAMPLING))
TAQ_USED_FILES=$(addsuffix _only.files,$(TAQ_USED_PATH))

# Create the necessary mothur .file files
$(TAQ_USED_FILES)\
$(PROC)/amp.files\
$(PROC)/mock_amp.files : code/make_amp_file.R
	R -e "source('code/make_amp_file.R')"


# Run the mothur analysis for error based on MOCK samples
$(PROC)/mock_error.% : code/mock_amp_mothur.batch
	bash code/mock_amp_mothur.batch


# Run the full analysis with different levels of subsampling with MOCK and Samples
$(PROC)/all_amp.% : code/amp_mothur.batch
	bash code/amp_mothur.batch

################################################################################
#
# Part 3: Metadata Processing and Analysis
#
#	Run scripts that analyze the generated data
#
################################################################################

# Create master meta data file
$(PROC)/amp.files : (TABLES)/meta_data.csv code/make_metadata_file.R
	R -e "source('code/make_metadata_file.R')"


# Set up fecal sample count files
FS_SHARED_PATH=$(addprefix $(PROC)/all_amp.0.03.subsample._,$(H_SAMPLING))
FS_SHARED_TABLES=$(addsuffix .shared,$(FS_PATH))
FS_PATH=$(addprefix $(TABLES)/fecal_sub_sample_,$(H_SAMPLING))
FS_COUNT_TABLES=$(addsuffix _count_table.csv,$(FS_PATH))
FS_ZSCORE_PATH=$(addprefix $(TABLES)/fecal_zscore_sub_sample_,$(H_SAMPLING))
FS_ZSCORE_TABLES=$(addsuffix _count_table.csv,$(FS_ZSCORE_PATH))


# Generate the fecal sample count tables
$(FS_SHARED_TABLES)\
$(TABLES)/meta_data.csv : $(FS_SHARED_TABLES) code/run_fecal_otu_count_tables.R
	R -e "source('code/run_fecal_otu_count_tables.R')"


# Generate the analysis tables for the number OTUs by subsampling for fecal samples
$(FS_COUNT_TABLES) : $(TABLES)/fecal_overall_anova_results.csv\
$(TABLES)/fecal_overall_tukey_results.csv $(FS_ZSCORE_TABLES)\
code/run_fecal_numOTU_analysis.R
	R -e "source('code/run_fecal_numOTU_analysis.R')"


# Set up mock sample count files
M_SHARED_PATH=$(addprefix $(PROC)/all_amp.0.03.subsample._,$(M_SAMPLING))
M_SHARED_TABLES=$(addsuffix .shared,$(M_SHARED_PATH))
M_PATH=$(addprefix $(TABLES)/mock_sub_sample_,$(M_SAMPLING))
M_COUNT_TABLES=$(addsuffix _count_table.csv,$(M_PATH))


# Generate the number of OTUs by subsampling for Mock samples
$(M_SHARED_TABLES) : $(TABLES)/mock_overall_anova_results.csv\
$(TABLES)/mock_overall_tukey_results.csv $(M_COUNT_TABLES)\
code/run_otu_diversity_analysis.R
	R -e "source('code/run_otu_diversity_analysis.R')"


# Set up mock seq error tables
M_ERROR_PATH=$(addprefix $(TABLES)/error_,$(M_SAMPLING))
M_ERROR_COUNT_TABLES=$(addsuffix _summary.csv,$(M_ERROR_PATH))
M_NUC_PATH=$(addprefix $(TABLES)/nucleotide_error_,$(M_SAMPLING))
M_NUC_TABLES=$(addsuffix _summary.csv,$(M_NUC_PATH))


# Generate the needed tables for error analysis and graphing
$(M_ERROR_COUNT_TABLES)\
$(M_NUC_TABLES) : $(PROC)/mock_error.count_table\
$(PROC)/mock_error.summary $(TABLES)/meta_data.csv\
code/run_seq_error_table_creations.R
	R -e "source('code/run_seq_error_table_creations.R')"

#Generate the analysis tables for the error metrics by subsampling for mocks
$(TABLES)/mock_error_overall_kruskal_results.csv\
$(TABLES)/mock_error_overall_dunn_results.csv\
$(TABLES)/mock_error_count_overall_kruskal_results.csv\
$(TABLES)/mock_error_count_overall_dunn_results.csv\
$(TABLES)/mock_chimera_overall_kruskal_results.csv\
$(TABLES)/mock_chimera_overall_tukey_results.csv : $(M_ERROR_COUNT_TABLES)\
code/run_error_analysis.R
	R -e "source('code/run_error_analysis.R')"

#Set up mock bray-curtis tables
M_BC_PATH=$(addprefix $(PROC)/all_amp.braycurtis.0.03.lt.,$(M_SAMPLING))
M_BC_TABLES=$(addsuffix .dist,$(M_BC_PATH))
M_BC_5_PATH=$(addprefix $(TABLES)/mock_bray_5_cycle_dist_,$(M_SAMPLING))
M_BC_5_TABLES=$(addsuffix _data.csv,$(M_BC_5_PATH))

#Set up human bray-curtis tables
FS_BC_PATH=$(addprefix $(PROC)/all_amp.braycurtis.0.03.lt.,$(H_SAMPLING))
FS_BC_TABLES=$(addsuffix .dist,$(FS_PATH))
FS_BC_5_PATH=$(addprefix $(TABLES)/bray_5_cycle_dist_,$(H_SAMPLING))
FS_BC_5_TABLES=$(addsuffix _data.csv,$(FS_BC_5_PATH))

#Generate the bray-curtis analysis for fecal samples
$(TABLES)/bray_sim_to_prev_cycle_kruskal_results.csv\
$(TABLES)/bray_permanova_by_taq_results.csv\
$(FS_BC_5_TABLES) : $(TABLES)/meta_data.csv $(FS_BC_TABLES)\
code/run_nmds_fecal_data.R
	 R -e "source('code/run_nmds_fecal_data.R')"

#Generate the bray-curtis analysis for stool samples
$(TABLES)/mock_bray_sim_to_prev_cycle_kruskal_results.csv\
$(TABLES)/mock_bray_permanova_by_taq_results.csv\
$(M_BC_5_TABLES) : $(TABLES)/meta_data.csv $(M_BC_TABLES)\
code/run_nmds_mock_data.R
	 R -e "source('code/run_nmds_mock_data.R')"


# Run the RF models for polymerase, cycles, and sample groups
$(PROC)/tables/sample_group_rf_imp_vars_summary.csv\
$(PROC)/tables/polymerase_group_top10_rf_imp_vars_summary.csv\
$(PROC)/tables/cycle_group_top10_rf_imp_vars_summary.csv\
$(PROC)/tables/polymerase_group_rf_model_summary.csv\
$(PROC)/tables/cycles_group_rf_model_summary.csv\
$(PROC)/tables/sample_group_rf_model_summary.csv\
$(PROC)/tables/polymerase_group_mock_rf_model_summary.csv\
$(PROC)/tables/cycles_group_mock_rf_model_summary.csv : code/run_rf_mock_test.R\
code/run_rf_test.R $(PROC)/all_amp.0.03.subsample.1000.shared $(PROC)/all_amp.taxonomy
	R -e "source('code/run_rf_mock_test.R')"
	R -e "source('code/run_rf_test.R')"


# Run the model testing
$(PROC)/tables/kruskal_rf_model_test.csv\
$(PROC)/tables/dunn_rf_model_test_summary.csv : $(PROC)/tables/sample_group_rf_model_summary.csv\
$(PROC)/tables/polymerase_group_rf_model_summary.csv\
$(PROC)/tables/cycles_group_rf_model_summary.csv\
$(PROC)/tables/polymerase_group_mock_rf_model_summary.csv\
$(PROC)/tables/cycles_group_mock_rf_model_summary.csv code/run_rf_model_probs_test.R
	R -e "source('code/run_rf_model_probs_test.R')"






################################################################################
#
# Part 4: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################

# GC content figure
$(FIGS)/Figure1.pdf : $(PROC)/gc_content_amp_summary.csv\
code/make_gc_content_graph.R
	R -e "source('code/make_gc_content_graph.R')"

# Run code to create Figure S2, 3, and 4 - Mock Sequences substitution frequency
$(FIGS)/FigureS1.pdf : $(PROC)/tables/full_nucleotide_error_summary.csv\
code/make_nucleotide_graphs.R
	R -e "source('code/make_nucleotide_graphs.R')"

# General Error graph
$(FIGS)/Figure2.pdf : $(M_ERROR_COUNT_TABLES) $(M_ERROR_COUNT_TABLES)\
$(M_ERROR_COUNT_TABLES) $(M_COUNT_TABLES) code/make_mock_chimera_count_graphs.R\
code/make_correlation_graphs.R
	R -e "source('code/make_mock_chimera_count_graphs.R')"

# Number of OTU graph
$(FIGS)/Figure3.pdf : $(FS_ZSCORE_TABLES) $(M_COUNT_TABLES)\
code/make_fecal_numOTU_graphs.R code/make_numOTU_graphs.R
	R -e "('code/make_numOTU_graphs.R')"


# Run code to create Figure 3 - Bray-Curtis distance differences by 5 set cycle
$(FIGS)/Figure4.pdf : $(M_BC_5_TABLES) $(FS_BC_5_TABLES)\
code/make_bray_distance_graphs.R
	R -e "source('code/make_bray_distance_graphs.R')"

# RF model classification for samples, cycles, and polymerase
$(FIGS)/Figure5.pdf : $(PROC)/tables/sample_group_rf_imp_vars_summary.csv\
$(PROC)/tables/polymerase_group_top10_rf_imp_vars_summary.csv\
$(PROC)/tables/cycle_group_top10_rf_imp_vars_summary.csv\
$(PROC)/sample_group_rf_model_summary.csv $(PROC)/polymerase_group_rf_model_summary.csv\
$(PROC)/cycles_group_rf_model_summary.csv $(PROC)/polymerase_group_mock_rf_model_summary.csv\
$(PROC)/cycles_group_mock_rf_model_summary.csv
	R -e "source('code/make_rf_fecal_imp_vars.graph.R')"



################################################################################
#
# Part 5: Pull it all together
#
# Render the manuscript
#
################################################################################


#write.paper : $(TABLES)/table_1.pdf $(TABLES)/table_2.pdf\ #customize to include
#				$(FIGS)/figure_1.pdf $(FIGS)/figure_2.pdf\	# appropriate tables and
#				$(FIGS)/figure_3.pdf $(FIGS)/figure_4.pdf\	# figures
#				$(FINAL)/study.Rmd $(FINAL)/study.md\
#				$(FINAL)/study.tex $(FINAL)/study.pdf
