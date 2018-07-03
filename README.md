## Assessing the Differences in 16S rRNA Gene Sequencing Due to High Fidelity DNA Polymerase 

**Background.** A typical 16S rRNA gene sequencing workflow can be divided into preservation, extraction, amplification, and sequencing steps. At each of these stages error can be introduced that will change the underlying bacterial community composition results. In this study we focus on the amplification step’s contribution to this overall error. To accomplish this we assessed 16S rRNA gene sequencing results in human fecal and mock community samples after using different high
fidelity (HiFi) DNA polymerases and number of amplification cycles.


**Methods.** We extracted DNA from fecal samples (n=4) using a PowerMag DNA extraction kit with a 10 minute bead beating step and amplified at 15, 20, 25, 30, and 35 cycles using Accuprime, Kappa, Phusion, Platinum, or Q5 HiFi DNA polymerase. Amplification of mock communities (technical replicates n=4) consisting of previously isolated whole genomes of 8 different bacteria used the same approach. We first assessed GC dependent differences, error rate, sequence error prevalence, chimera prevalence, and correlation between chimera prevalence and number of Operational Taxonomic Units (OTUs) by polymerase and number of cycles. Next, differences in the number of OTUs and taxa was examined based on the polymerase and number of cycles used. Additionally, differences in the bacterial community composition by the Bray-Curtis index also was assessed based on polymerase and number of cycles. Finally, Random Forest models were created to test whether the bacterial community was better at classifying polymerases, number of cycles, or
individual donor. We also assessed whether the most important taxa in the polymerase and number of cycle Random Forest models also were the most important in the model for individual donors.


**Results.** Although the total counts for specific taxa varied by polymerase, we found no noticeable differences in total counts based on high and low GC content (P-value $\geqslant$ 0.23). Chimera prevalence in mock communities varied by polymerase with differences being most notable at 35 cycles (Kappa = 5.71% (median) versus Platinum = 26.62%) and this variation persisted after chimera removal using VSEARCH. We also observed positive correlations between chimera prevalence and the number of OTUs with Platinum having the highest (R<sup>2</sup> = 0.974) and Kappa having the worst (R<sup>2</sup> = 0.478). When analyzing mock community samples the variation in the number of OTUs detected by the polymerases was observable as early as 20 cycles (P-value = 0.002). There also was a large range in the number of OTUs amplified by the polymerases at 35 cycles (Accuprime = 15 - 20 versus Phusion = 14 - 73). When analyzing fecal samples we observed that the range in the number of OTUs detected was not consistent between HiFi DNA polymerases (eg. at 35 cycles Accuprime = 84 - 106 (min - max) versus Phusion = 84 - 136). Additionally, the median number of OTUs varied by HiFi DNA polymerase used at 35 cycles (P-value < 0.0001). Random Forest models were most successful at classifying individual donor samples rather than polymerase or number of cycles used (P-value $\leqslant$ 5.49e-07). Additionally, the most important OTUs in the polymerase and number of cycle models were not the most important in the individual donor sample model.


**Conclusions.** Although there are 16S rRNA gene sequencing differences based on polymerase and number of cycles used, they are small with respect to the biological differences between indivdiuals. Additionally, our observations were not able to identify GC-based differences in amplification based on the number of cycles used. Collectively, these results provide evidence that a real biological difference between groups, based on 16S rRNA gene sequencing, should be detectable regardless of polymerase and number of cycles used.





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
