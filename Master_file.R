mychoice <- menu( c("No","Yes"), graphics=TRUE, title="Install required libraries?" )

if (mychoice==2){ 
	install.packages("MASS")
	install.packages("RcppNumerical")
	install.packages("BMS")
	install.packages("BAS")
	install.packages("parallel")
	install.packages("fda.usc")
	install.packages("coda")
} 

source("./Files/MAdaSub_LOAD.R")

choices <-c( "Illustrative simulated example (Section 5.1) ", 
             "Simulation study (Section 5.2) ", 
             "Tecator data application (Section 6.1) ", 
             "PCR data application (Section 6.2) ", 
             "Leukemia data application (Section 6.3) ")
mychoiceM <- menu( choices , graphics=TRUE, title="Which setting \n do you want to run?" )

if (mychoiceM==1) source("./Files/MAdaSub_Low_gprior_with_correlation.R")
if (mychoiceM==2) source("./Files/MAdaSub_Simulations.R")
if (mychoiceM==3) {
  choices <-c( "Multiple serial and parallel MAdaSub chains ",
               "One serial MAdaSub chain ")
  mychoiceMa <- menu( choices , graphics=TRUE, title="Which MAdaSub chains should be used?" )
  
  if (mychoiceMa==1) source("./Files/MAdaSub_parallel_Tecator.R")
  if (mychoiceMa==2) source("./Files/MAdaSub_Tecator_data.R")
}

if (mychoiceM==4){  
	# Example on PCR data 
	# Please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)
	# Download file from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip
	# extract  files Xgenes.txt, gene_id.txt & Y3.txt from 77-5.Song.zip inside the zip file
	print( noquote("STEP 1: Please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)" ))
	print( noquote("        from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip" ))
	print( noquote("STEP 2: Extract  files Xgenes.txt, gene_id.txt & Y3.txt from 77-5.Song.zip from the zip file to the subfolder 'Files' "))
	print( noquote("        (path: '77-5.Song.zip/data_code.zip/data_and_code/file/data/PCR')" ))

	menu("Press 1 when you are ready")
	
	choices <-c( "Multiple serial and parallel MAdaSub chains ",
	             "Three serial MAdaSub chains ")
	mychoiceMa <- menu( choices , graphics=TRUE, title="Which MAdaSub chains should be used?" )
	
	if (mychoiceMa==1) source("./Files/MAdaSub_parallel_PCR.R")
	if (mychoiceMa==2) source("./Files/MAdaSub_PCR_data.R")
	
} 

if (mychoiceM==5) {
  mychoice <- menu( c("No","Yes"), graphics=TRUE, title="Install required libraries?" )
  
  if (mychoice==2){ 
    install.packages("BiocManager")
    BiocManager::install("genefilter")
    BiocManager::install("golubEsets")
  } 
  
  choices <-c( "Multiple serial and parallel MAdaSub chains ",
               "Three serial MAdaSub chains ")
  mychoiceMa <- menu( choices , graphics=TRUE, title="Which MAdaSub chains should be used?" )
  
  if (mychoiceMa==1) source("./Files/MAdaSub_parallel_Golub.R")
  if (mychoiceMa==2) source("./Files/MAdaSub_Golub_data.R")
}





