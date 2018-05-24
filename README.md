## The Influence of High Fidelity DNA Polymerase on 16S rRNA Gene Sequencing

**Background.** Bias is introduced at different stages of the 16S rRNA gene sequencing workflow. Although the number of cycles and high fidelity (HiFi) DNA polymerase are studied less often, they are still important sources of bias in this workflow. Here, we examine how both cycle number and polymerase can change the bacterial community and introduce both error and bias to the final obtained results.


**Methods.** We extracted DNA from fecal samples (n=4) using a PowerMag DNA extraction kit with a 10 minute bead beating step and amplified at 15, 20, 25, 30, and 35 cycles using Accuprime, Kappa, Phusion, Platinum, or Q5 HiFi DNA polymerase. Amplification of mock communities (technical replicates n=4) consisting of previously isolated whole genomes of 8 different bacteria used the same approach. The analysis initially examined the number of Operational Taxonomic Units (OTUs) for fecal samples and mock communities. It also assessed polymerase dependent differences in the Bray-Curtis index, error rate, sequence error prevalence, chimera prevalence, and chimera prevalence correlation with the number of OTUs.



**Results.** When analyzing fecal samples we observed a different number of OTUs between HiFi
DNA polymerases at 35 cycles (P-value < 0.0001). Our analysis identified these polymerase
dependent differences in the number of OTUs as early as 20 cycles in mock communities (P-value
= 0.002). Chimera prevalence varied by polymerase and this variation persisted after chimera
removal using VSEARCH. We also observed positive correlations between chimera prevalence
and the number of OTUs.



**Conclusions.** HiFi DNA polymerase dependent differences in the number of OTUs and chimera
prevalence makes comparison across studies difficult. Care should be exercised when choosing
the polymerase and number of cycles to be used in 16S rRNA gene sequencing studies.




### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make (v3.81) should be located in the user's PATH
* mothur (v1.39.3) should be located in the user's PATH
* sratoolkit (2.8.2-1) should be located in the user's PATH
* R (v.3.4.2) should be located in the user's PATH


#### Running analysis

```
git clone https://github.com/SchlossLab/Sze_PCRSeqEffects_XXXX_2017.git
make write.paper
```
