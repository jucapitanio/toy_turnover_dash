# Code example for Dr. Sagar Damle: toy_turnover_dash
Simplified example code with toy data for Dr. Sagar Damle to evaluate before interview.

### What it does:
It analyzes TMT-SILAC data to estimate the rate of synthesis and degradation of proteins/peptides using a couple different methods of non-linear regressions.

### How to use:
The easiest way to use this is to open the dashboard file dashboad_mockup.Rmd in RStudio, edit the 2 required experimental parameters at the top (example parameters provided) and knit the markdown file.
This will perform all the analysis, create a dashboard to visualize the results, and save the outputs in a new folder.
An option for running only the script behind the dashboard in the command line or as a sourced local job in RStudio is also described in the script itself (toy_data_analysis.R).

### Inputs:
The input consists of 3 tables, SILAC quantification data (MS1), TMT quantification of light labelled peptides and TMT quantification of heavy labelled peptides (both MS3 and NOT normalized). These quantifications were originally done using the IP2 proteomics pipeline. See toy_data/IP2_data for example datasets. I've added 2 examples, one from the nuclear and one from the cytoplasmic proteome. Due to confidentiality these have been filtered to contain only peptides from the 60S ribosomal proteins, so results are not very representative. 

### Outputs
The output is an html dashboard with multiple tabs and subtabs that can be viewed in any browser (dashboad_mockup.html). 
Additional tables of the processed data and model results are saved in toy_data/output_dir. The results are also saved as R objects in toy_data/scripts_results. 

#### Notes: 
- As I mentioned, this is a simplified older version of the project. You'll notice that in this older commit part of the functions are properly documented with roxygen headers and part are not since I was in the process of building the R package. In the final version we have additional models, error estimates, and model tests. In the new version all the functions are also built into an R package that is distributed (along with the other package dependencies) in a docker image. 
- Please do not share or utilize any of the code or data present here without explicit written authorization by Dr. Juliana Capitanio and Prof. Martin Hetzer

Thank you and please let me know if you have any questions.

Juliana