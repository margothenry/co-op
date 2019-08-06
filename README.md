Recommended workflow:

_ Learn what a suitable paper for this analysis contains: Go to "docs/Finding a paper". Read the newest version. 
_ Once you've read the paper & determined that it's useful, download the paper (in PDF) to papers_toAnalyze. Naming convention is 1st author's last name followed by year e.g. Tenaillon2016.
_ Download the data to "data_in/to_analyze", keeping the data in a folder named identically with its paper (AuthorYear).
_ The data file as downloaded should be named "AuthorYear_original".
_ Do some data cleaning in Excel. Recommended data cleaning steps can be found in "docs/Recommended data cleaning steps."
_ Once done, save-as the data in .csv format ("CSV (Comma delimited) (*.csv)") & name the new file "AuthorYear_usable".
_ Make sure there are no space between columns in the _usable.csv file.
_ Move the folder containing the data from "to_analyze" to "original & usable".
_ Open "R/Master_Data". Read the existing codes. Pay extra attention to the first few datasets of each function - those are usually notated at every step to explain its purpose.
_ Analyze the data. Depending on the data format, use the appropriate function call. If unclear about which function call to use, open "R/dgconstraint/functions" & read each R file.
_ In papers_done, create a folder named with your initials (first & last names).
_ Add an underscore & your initials to the name of the paper you just analyzed e.g. Tenaillon2016_MH. 
_ Move the paper from papers_toAnalyze to your sub-folder in papers_done. 
_ Commit to Github after every change made to the "co-op" folder.
_ At the end of the week, push the "co-op" folder to Github.
_________________________________________

A run-down of each folder in "co-op":

_ data_in: Input for functions, in various stages of cleaning: to_analyze -> original & usable -> for_func.
--- for_func: The data that results from cleaning "original & usable" in R.
_ data_out: Output of functions.
--- analyses: Contains individual .csv files of each dataset analyzed.
--- images: Some screenshots & visualizations of the output. 
--- intermediate: By-products of the functions. You don't need to worry about this.
--- master_analyses: A .csv file containing analyses of all datasets.
_ docs: Some documents to help you with your analyses.
--- Finding a paper: Contains different versions of the eponymous document, written to help you look for suitable papers for the analysis.
--- Paper errors: Mistakes spotted in the papers you have read. 
--- Recommended data cleaning steps: Helps you clean data in Excel & R.
_ paper_screening: Info about paper screening. As of 190806, the folder includes an overview of the papers screened.
_ papers_done: All papers analyzed. Subfolders are named after the first & last name initials of the analyst.
_ papers_toAnalyze: Papers deemed useful, but whose data has yet to be analyzed.
_ R: Contains R scripts.
--- dgconstraint: 
------ functions: Scripts for functions used to analyze different kinds of datasets.
------ inst: Contains "GeneDatabase", a file with the different species included in the analyses & the corresponding number of genes.
------ vignette: Contains a detailed Rmarkdown guide to the different functions in dgconstraint, as well as the figures in the guide (folder "figs").
--- Data_Visualization: The codes used to generate the visualizations in "data_out/images".
--- Master_Data: The codes used to analyze the data.
--- StatisticalAnalysis: Some quantitative analysis that might be helpful.