# Steps to This Analysis

## 1. Ensure that BioLockJ is on local system
BioLockJ - https://biolockj-dev-team.github.io/BioLockJ/Getting-Started/

## 2. Download the PublicationAnalysis directory
git clone https://github.com/FarnazFouladi/Anorexia_Analyses2021

## 3. Set up R and required packages

Make sure R is installed.  See https://www.r-project.org/.  These scripts were written with 4.0.2.

Make sure all required R packages are installed (see _load.packages.R_).  

To do this, make sure the _load.packages.R_ script runs without errors.

Move to the analysis folder:            
`cd Anorexia_Analyses2021`
`cd analysis/`

Run the library list script.
`Rscript resources/load.packages.R `

## 4. Run BioLockJ pipeline

Move to the analysis folder:            
`cd Anorexia_Analyses2021`
`cd analysis/`

To run the pipeline using **locally installed software**:                 
`biolockj AnorexiaAnalysis.config`
