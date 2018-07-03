## Assessing differences in high fidelity DNA polymerase for 16S rRNA gene sequencing

**Background.** It is challenging to compare 16S rRNA gene sequencing data across studies and
one of the reasons for this is due to error. There are many different places throughout the workflow where error can be introduced into the pipeline. Here, we focus on studying how the number of cycles and high fidelity (HiFi) DNA polymerase introduce error by varying cycle number and polymerase used to amplify 16S rRNA genes in human fecal and mock community samples.



**Methods.** We extracted DNA from fecal samples (n=4) using a PowerMag DNA extraction kit with a 10 minute bead beating step and amplified at 15, 20, 25, 30, and 35 cycles using Accuprime, Kappa, Phusion, Platinum, or Q5 HiFi DNA polymerase. Amplification of mock communities (technical replicates n=4) consisting of previously isolated whole genomes of 8 different bacteria used the same approach. The analysis initially examined the number of Operational Taxonomic Units (OTUs) for fecal samples and mock communities. It also assessed polymerase and number of cylce dependent differences in the Bray-Curtis index, error rate, sequence error prevalence, chimera prevalence, and the correlation between chimera prevalence and number of OTUs.



**Results.** When analyzing fecal samples we observed that the range in the number of OTUs
detected was not consistent between HiFi DNA polymerases (eg. at 35 cycles Accuprime = 84 - 106
(min - max) versus Phusion = 84 - 136). Additionally, the median number of OTUs varied by HiFi
DNA polymerase used 35 cycles (P-value < 0.0001). When analyzing mock community samples the variation in the number of OTUs detected by the polymerases was observable as early as 20 cycles (P-value = 0.002). There also was a large range in the number of OTUs amplified by the polymerases at 35 cycles (Accuprime = 15 - 20 versus Phusion = 14 - 73). Chimera prevalence in mock communities varied by polymerase with differences being most notable at 35 cycles (Kappa = 5.71% (median) versus Platinum = 26.62%) and this variation persisted after chimera removal using VSEARCH. We also observed positive correlations between chimera prevalence and the number of OTUs with Platinum having the highest (R<sup>2</sup> = 0.974) and Kappa having the worst (R<sup>2</sup> = 0.478).




**Conclusions.** Although the variation in the number of OTUs in fecal samples could be due to
certain polymerases capturing the biological variability better than others, this is unlikely to be the main reason for our observed differences. In mock community samples, the strong correlation between chimera prevalence and the number of OTUs suggests that this is the main reason for differences between the polymerases. Ultimately, this variation makes comparison across studies difficult and care should be exercised when choosing the polymerase and number of cycles in 16S rRNA gene sequencing studies.





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
