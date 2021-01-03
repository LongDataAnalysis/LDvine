# LDvine

This repository contains the scripts, data, results, and examples for the project "A time-heterogeneous D-vine copula model for unbalanced and unequally spaced longitudinal data".

The repository is structured as follows:

-`data`: Contains csv/rda/rds/txt files that are used as data into R in order to do the analysis. 

	We used a subset of Manitoba Follow-up Study (MFUS) data. Because of data confidentiality, the data set file is not included in this R project. 
	However, the data descriptions and structures are presented below:   

	VARIABLE NAME AND DESCRIPTIONS:

	COL 01 -- id = individual id number

	COL 02 -- followup = follow-up number of each individual

	COL 03 -- px_age = age at the time of experiment

	COL 04 -- bmi = body mass index in each follow-up

	COL 05 -- IHD = ischemic heart disease status

	COL 06 -- DBP = Diastolic Blood Pressure

	Finally, the data set has 462 subjects with at most 5 follow-ups. In total, the data has 1713 rows and 6 columns.

- `src`: Contains R functions (.R) that implement the methods in this R project. 

- `examples`: Contains RMarkdown (.Rmd) files with the complied pdf versions. This folder contains two sub-folders:	

	- `Simulation`: this sub-folder contains .Rmd files with the complied pdf versions for simulation study.
RMarkdown files: The R-code for the simulation example is in the R Markdown file "LD-vine_Sim.Rmd" 
(use RStudio for opening/modifying/compiling). The compiled version is available as a pdf ("LD-vine_Sim.pdf"). 
The code is commented and directly evaluated (see "LD-vine_Sim.pdf"). 
To reproduce the simulated example,  simply run the R code line by line (after having installed all R-packages 
that are loaded at the very beginning - if needed).
	- `Data Analysis`: this sub-folder contains the data application code in RMarkdown. Only the complied version is available as a pdf file ("LD-vine_Data.pdf"), which includes the partial results of data analysis presented in Section 5 of manuscript. However, since the real data is not provided because of data confidentiality, this RMarkdown file is not reproducible.   

- `results`: Save all printed project outputs here, including plots, HTML, and data exports.


## First time only: Setting up the directory

To start you need to open the project file (LDvine.Rproj) in RStudio which will bring RStudio to the root directory.




