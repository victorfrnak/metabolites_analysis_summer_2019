install.packages("VGAM")
library(VGAM)
#SorghumWTF

#fileSorghum = read.csv(file = "/Users/ahubbard/sorghum_least_stringent_settings_first_analysis.txt", sep = "\t", header = TRUE, comment.char = "",quote = "",stringsAsFactors = FALSE, na.strings = "NA")


#setwd("/Users/ahubbard/xcms_metabolomics/working_with_most_realistic_data")
#setwd("/Users/ahubbard/xcms_metabolomics/wrapping_up_for_ivan_sorghum")

#setwd("/Users/ahubbard/xcms_metabolomics/wrapping_up_for_ivan_sorghum_retry_with_may10th")
#setwd("/Users/ahubbard/xcms_metabolomics/wrapping_up_for_ivan_sorghum_retry_with_may12th")
setwd("/Users/ahubbard/xcms_metabolomics/sorghum_final_may13th")

fileSorghum = read.csv(file = "/Users/ahubbard/sorghum_more_stringent_settings_minfrac_and_ppm_first_analysis_to_read_in.txt", sep = "\t", header = TRUE, comment.char = "",quote = "",stringsAsFactors = FALSE, na.strings = "NA")
namesForFileSorghum = paste(fileSorghum$mz, fileSorghum$rt, sep = "_")
rownames(fileSorghum) = namesForFileSorghum
#rownames(fileSorghum) = make.names(fileSorghum$Metabolite.Name, unique = TRUE)

fileSorghum$X = NULL
fileSorghum$mz = NULL
fileSorghum$rt = NULL


#change the names if needed
for(i in 1:length(colnames(fileSorghum)))
{
  
  
  colnames(fileSorghum)[i] = gsub(".y", "",as.character(colnames(fileSorghum)[i]))
  colnames(fileSorghum)[i] = gsub(".x", "",as.character(colnames(fileSorghum)[i]))
  #gsub('.x', colnames(fileSorghum)[i])
  toTest = strsplit(as.character(colnames(fileSorghum)[i]), "_")[[1]][length(strsplit(as.character(colnames(fileSorghum)[i]), "_")[[1]])]
  
  #print(i)
  #print(toTest)
  
  
  #if we need to, remove the last field
  if(grepl("T", as.character(toTest)) == FALSE)
  {
    toReplace = paste(strsplit(as.character(colnames(fileSorghum)[i]), "_")[[1]][-length(strsplit(as.character(colnames(fileSorghum)[i]), "_")[[1]])], collapse  = "_")
    print(toReplace)
    
    colnames(fileSorghum)[i] = toReplace
    
  }
  #remove the .x and .y
  
}


#now, we'll join the transcriptome data
#read in the the transcriptome metadata
transcriptMetaSorghum = read.table(file = "/Users/ahubbard/rna_seq_metadata_sorghum_preliminary.txt",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
transcriptMetaSorghum$Sample = gsub("(.*_P1_)","",transcriptMetaSorghum$Sample)
transcriptMetaSorghum$Sample = gsub("(.*_P2_)","",transcriptMetaSorghum$Sample)


### read in the metadata ###
designSetupSorghum = read.table("/Users/ahubbard/Downloads/sorghum_sample_list_to_read_in.txt", sep = "\t", header = TRUE)
designSetupSorghumWithTubeinfo = read.table("/Users/ahubbard/Downloads/sorghum_sample_list_to_read_in.txt", sep = "\t", header = TRUE)

#name of metabolite samples with RNA seq data as well
metaboliteSamplesWithRNASeq = designSetupSorghumWithTubeinfo[designSetupSorghumWithTubeinfo$Tube[intersect(as.numeric(designSetupSorghumWithTubeinfo$Tube), as.numeric(transcriptMetaSorghum$Sample))], ]$Name


#store all of the timepoints data in a vecotr
timePointsToParseSorghum = as.vector(unique(designSetupSorghum$Time))

#don't cut out T0 this time
#timePointsToParseSorghum = timePointsToParseSorghum[-1]




#vec of libs we'll want to pull out
#this is going to be metabolite data for which we
#also have transcriptome info
#forPCATPMHeadersSorghum = read.table("/Users/ahubbard/drought_stress_joined_and_sorted_sorghum_FPKM_tables_large_to_cluster_should_have_all.txt", sep = "\t", header = TRUE)
forPCATPMHeadersSorghum = read.table("/Users/ahubbard/xcms_metabolomics/to_test_with_rejoin.txt", sep = "\t", header = TRUE)


colnames(forPCATPMHeadersSorghum) = gsub(".*_P1_","",colnames(forPCATPMHeadersSorghum))
colnames(forPCATPMHeadersSorghum) = gsub(".*_P2_","",colnames(forPCATPMHeadersSorghum)) 

rownames(forPCATPMHeadersSorghum) = forPCATPMHeadersSorghum$Name
forPCATPMHeadersSorghum$Name = NULL

#store justthe number of the tube from the RNA seq metadata
colnames(forPCATPMHeadersSorghum) = gsub("(.*)_.*_.*_.*_.*", "\\1", colnames(forPCATPMHeadersSorghum))
colnames(forPCATPMHeadersSorghum) = as.numeric(colnames(forPCATPMHeadersSorghum))


#convert RNA seq data to numeric 
forPCATPMHeadersSorghum[] <- lapply(forPCATPMHeadersSorghum, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})


#Now loop through the timepoint and do the analysis from start
#to finish
for(i in 1:length(colnames(forPCATPMHeadersSorghum)))
{
  tubeNumber = colnames(forPCATPMHeadersSorghum)[i]
  
  
  #look up the full sample name of each RNA seq sample
  #in the metadata table
  vecOfPossibleNames = designSetupSorghumWithTubeinfo[designSetupSorghumWithTubeinfo$Tube == tubeNumber,]
  nameToAdd = as.character(vecOfPossibleNames$Name)
  
  print(nameToAdd)
  
  colnames(forPCATPMHeadersSorghum)[i] = nameToAdd
  
}




#if you want to subset ....
#forPCATPMHeadersSorghum = forPCATPMHeadersSorghum[rownames(forPCATPMHeadersSorghum) %in% intersect(noquote(rownames(forPCATPMHeadersSorghum)), as.vector(genesWithArabidopsisHits[,1])),]
#we will now have consistent naming with the metabolite data
#so that we can find the corresponding transcriptome data

vectorOfTranscriptInfo = vector()
RNAlibraryNumbers = vector()


theCombinedGenesAndMetabNames = c(rownames(fileSorghum),rownames(forPCATPMHeadersSorghum))
infoOnEachMetaboliteSorghum <- data.frame(zeros(matrix(ncol = 2, nrow = length(theCombinedGenesAndMetabNames))))
colnames(infoOnEachMetaboliteSorghum) = c("Comparisons_drought_responsive","Which Time Points")
rownames(infoOnEachMetaboliteSorghum) = theCombinedGenesAndMetabNames



#have a table to store the drought responsive metabolites
#and drought responsive 


parentDir = getwd()
responsiveInfoTableSorghum <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParseSorghum)),y=0)) 
colnames(responsiveInfoTableSorghum) = c("responsive_metabolites","responsive_genes")
rownames(responsiveInfoTableSorghum) = timePointsToParseSorghum


plsDAInfoTable <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParseSorghum)),y=0)) 
colnames(plsDAInfoTable) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTable) = timePointsToParseSorghum

#store the pls-da for just metabs
plsDAInfoTablejustMetsSorghum <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParseSorghum)),y=0)) 
colnames(plsDAInfoTablejustMetsSorghum) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTablejustMetsSorghum) = timePointsToParseSorghum

plsDAInfoTableJustRNA <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParseSorghum)),y=0)) 
colnames(plsDAInfoTableJustRNA) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTableJustRNA) = timePointsToParseSorghum



#Metabolites, classifying by all three gentoypes and drought treatment
plsDAInfoTableJustMetabAllThree <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParseSorghum)),y=0)) 
colnames(plsDAInfoTableJustMetabAllThree) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTableJustMetabAllThree) = timePointsToParseSorghum


#have a table keeping track of all the plsda permutation info
plsDAInfoTablePermutations <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParseSorghum)),y=0)) 
colnames(plsDAInfoTablePermutations) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTablePermutations) = timePointsToParseSorghum





#put everything in here, this directory
#setwd("/Users/ahubbard/xcms_metabolomics")


fileSorghum[fileSorghum == 0] <- NA
fileSorghum = na.omit(fileSorghum)

tableNoMissingSorghum = fileSorghum


for( t in 1:(length(timePointsToParseSorghum)))
{
  
  #just make sure that we re-set the table each time
  tableNoMissingSorghum = fileSorghum
  fileToJoin = t(tableNoMissingSorghum)
  
  timeSample = as.vector(unique(timePointsToParseSorghum)[t])
  
  pdf(paste(timeSample, "allPlotsMay12th.pdf", sep = "_"))
  
  #subset by just the time
  theSubsetSorghum = subset(designSetupSorghum, Time == timeSample)
  
  theDirForOutput = paste0(parentDir, "/",paste(timeSample,collapse = ""),"_full_report_all_one_file")
  
  #store everything in a single output directory 
  mkdir(theDirForOutput)
  
  #because we'll take out the special characters later
  #but need to match with genotype now
  mkdir(paste0("/Users/ahubbard/PCA_testing_all_Sorghum_no_outliers/", timeSampleForDir))
 
  
  #get a subset od all the libraries 
  #in this sample
  theSubsetSorghum = subset(designSetupSorghum, Time %in% timeSample)
  tableNoMissingSorghum = tableNoMissingSorghum[,colnames(tableNoMissingSorghum) %in% theSubsetSorghum$Name]
  tableNoMissingSorghum = tableNoMissingSorghum
  
  tableNoMissingSorghum = as.data.frame(tableNoMissingSorghum)
  
  justMetabsSorghum = tableNoMissingSorghum
  
  
  #seems that we have factors where we want numbers
  #but we don't want the NA's to be zeroes
  tableNoMissingSorghum[] <- lapply(tableNoMissingSorghum, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  

  #make sure that we have a copy of the tableNoMissingSorghum 
  #for just the metabolites in order to do the plsa-da
  
  #justMetabsSorghum = tableNoMissingSorghum
  

  #focus only on those for which we have corresponding metabolite data
  tableNoMissingSorghum = tableNoMissingSorghum[,colnames(tableNoMissingSorghum) %in% colnames(forPCATPMHeadersSorghum)]
  toMergeWithMetabsSorghum = forPCATPMHeadersSorghum[,colnames(forPCATPMHeadersSorghum) %in% colnames(tableNoMissingSorghum)]
  
  justMetabsSorghum = tableNoMissingSorghum
  

  tableNoMissingSorghum = rbind(tableNoMissingSorghum, toMergeWithMetabsSorghum)
  

  accessionPCA = paste0(theDirForOutput, "/allTimepoinByAccession.jpg")
 
  #now we will combine the metabs and transcriptome data
  ##toSubsetFromTranscriptome = colnames(tableNoMissingSorghum)
  
  #focus only on those for which we have corresponding metabolite data
  
  #tableNoMissingSorghum = tableNoMissingSorghum[,colnames(tableNoMissingSorghum) %in% colnames(forPCATPMHeadersSorghum)]
  #toMergeWithMetabsSorghum = forPCATPMHeadersSorghum[,colnames(forPCATPMHeadersSorghum) %in% colnames(tableNoMissingSorghum)]
  
  
  
  #log transform the table of metabolites and rbind with gene expression
  tableNoMissingSorghum = log(tableNoMissingSorghum + .001)
  tableNoMissingSorghum = rbind(toMergeWithMetabsSorghum, tableNoMissingSorghum)
  
  
  forPCANoZeroes = tableNoMissingSorghum
  #rownames(forPCANoZeroes) = dataAndMetadataMerged$Name
  
  #before we do the tests, we want to convert to numeric 
  #run the pls-da and keep track
  #of the outputs
  
  
  treatmentPCA = paste0(theDirForOutput, "/allTimepointByTreatment.pdf")
  
  #also, do the metabolome and the transcriptome separately
  XjustMetsSorghum = t(justMetabsSorghum)
  YjustMetsSorghum = rownames(t(justMetabsSorghum))
  
  YjustMetsSorghum = gsub(".*PI_329311","",YjustMetsSorghum)
  YjustMetsSorghum = gsub(".*Grassl","",YjustMetsSorghum)
  YjustMetsSorghum = gsub(".*BT_623","",YjustMetsSorghum)
  YjustMetsSorghum = gsub("Rep_.*","",YjustMetsSorghum)
  YjustMetsSorghum = substr(YjustMetsSorghum,1,nchar(as.character(YjustMetsSorghum))-1)
  
  
  set.seed(2543) 
  plsdatestjustMetsSorghum <- plsda(as.matrix(XjustMetsSorghum), as.factor(YjustMetsSorghum), ncomp = 2)
  plotIndiv(plsdatestjustMetsSorghum, ind.names = TRUE, ellipse = TRUE, legend = TRUE, mainTitle = paste0("PLS-DA at Combined Metabolome + Transcriptome", timePointsToParseSorghum))
  
  
  set.seed(2543) 
  perf.plsdajustMetsSorghum <- perf(plsdatestjustMetsSorghum, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = 10)
  #plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal", mainTitle = paste0("PLS-DA at Just Metabolome", timePointsToParseSorghum))
  
  treatmentPCA = paste0(theDirForOutput, "/allTimepointByTreatment.pdf")

  myPLSDAErrorJustMetInfo =  perf.plsdajustMetsSorghum$error.rate$overall[rownames(perf.plsda$error.rate$overall) == "comp 2",colnames(perf.plsda$error.rate$overall) == "centroids.dist"]
  
  #store the output of myPLSDAError, note well !!! 
  plsDAInfoTablejustMetsSorghum[rownames(plsDAInfoTablejustMetsSorghum) == timeSample, colnames(plsDAInfoTablejustMetsSorghum) == "plsdaInfo"] = myPLSDAErrorJustMetInfo
  plsDAInfoTablejustMetsSorghum[rownames(plsDAInfoTablejustMetsSorghum) == timeSample, colnames(plsDAInfoTablejustMetsSorghum) == "bestCompNumber"] = which.min(perf.plsda$error.rate$overall)
  
  
  #now do just the transcriptome
  XJustTransSorghum = t(forPCATPMHeadersSorghum)
  YJustTransSorghum = rownames(t(forPCATPMHeadersSorghum))
  
  YJustTransSorghum = gsub(".*PI_329311","",YJustTransSorghum)
  YJustTransSorghum = gsub(".*Grassl","",YJustTransSorghum)
  YJustTransSorghum = gsub(".*BT_6231","",YJustTransSorghum)
  #YJustTransSorghum = gsub("Rep_.*","",Y)
  
  
  YJustTransSorghum = substr(YJustTransSorghum,1,nchar(as.character(YJustTransSorghum))-1)
  
  
  set.seed(2543) 
  plsdatestSorghumJustTrans <- plsda(as.matrix(XJustTransSorghum), as.factor(YJustTransSorghum), ncomp = 2)
  plotIndiv(plsdatestSorghumJustTrans, ind.names = TRUE, ellipse = TRUE, legend = TRUE, mainTitle = paste0("PLS-DA at Combined Metabolome + Transcriptome", timePointsToParseSorghum))
  
  
  #focus on only the compounds that are drought responsive
  #at one or more time points
  topFeatureInfoSorghum = vector()  
  
  #pooling the genotypes together into just control an drought stress
  tableNoMissingSorghum = justMetabsSorghum
  for (a in 1:nrow(tableNoMissingSorghum))
  {
    
    theCompound = rownames(tableNoMissingSorghum)[a]
    
    #  timePoint = timeSample 
    vec40Percent = vector()
    vec100Percent = vector()
    
    #theSubsetSorghum = subset(designSetupSorghum, Time %in% timeSample)
    for (k in 1: (length(unique(designSetupSorghum$Treatment))))
    {
      
      treatmentSample = unique(designSetupSorghum$Treatment)[k]  
      
      theSubsetSorghumII = theSubsetSorghum[theSubsetSorghum$Treatment == treatmentSample,]
      theSamplesToImputeSorghum = theSubsetSorghumII$Name
      
      ### get the info for the metabolite ###
      forTtest = tableNoMissingSorghum[rownames(tableNoMissingSorghum) == theCompound,intersect(colnames(tableNoMissingSorghum), theSamplesToImputeSorghum)]
      
      
      if(treatmentSample == 40)
      {
        vec40Percent = forTtest
      }
      
      if(treatmentSample == 100)
      {
        vec100Percent = forTtest
      }
    }
    
    if(length(unique(as.numeric(vec40Percent)) > 1) && (length(unique(as.numeric(vec100Percent))) > 1))
    {
      theTest = t.test(as.numeric(vec40Percent), as.numeric(vec100Percent))
      pvalue = theTest$p.value
    }
    else
    {
      pvalue = 1
    }
    
    if(pvalue < .05)
    {
      print("we have a hit !!!")
      infoOnEachMetaboliteSorghum[rownames(infoOnEachMetaboliteSorghum) == theCompound, 1] =  infoOnEachMetaboliteSorghum[rownames(infoOnEachMetaboliteSorghum) == theCompound, 1] + 1
      
      #tTestResults = theTest$estimate
      if(tTestResults[1] > tTestResults[2])
      {
        resultsForTable = paste0("upInDrought",timeSample) 
      }
      if(tTestResults[2] > tTestResults[1])
      {
        resultsForTable = paste0("downInDrought",timeSample) 
      }
      
      infoOnEachMetaboliteSorghum[rownames(infoOnEachMetaboliteSorghum) == theCompound, 2] =  paste(infoOnEachMetaboliteSorghum[rownames(infoOnEachMetaboliteSorghum) == theCompound, 2], resultsForTable, collapse = ",")
      topFeatureInfoSorghum = c(topFeatureInfoSorghum, theCompound)
      
      #we'll also want to do by genotype, I assume 
    }
  }
  
  fileForTopFeatures = paste0(parentDir,"/drought_responsive_metabolites_and_genes_may12th", paste(timeSample,collapse = ''), ".txt", sep ="_")
 
  #output just the genes
  fileForTopFeaturesjustMetabsSorghumSorghum = paste0(parentDir,"/drought_responsive_metabolites_only", paste(timeSample,collapse = '_'), ".txt")
  topFeatureInfoSorghumjustMetabsSorghum = topFeatureInfoSorghum[topFeatureInfoSorghum %in% rownames(tableNoMissingSorghum)]
  write.table(topFeatureInfoSorghumjustMetabsSorghum, file = fileForTopFeaturesjustMetabsSorghumSorghum, row.names =  FALSE)
  
  topMetabRounded = round(as.numeric(gsub("(.*)_.*", "\\1", topFeatureInfoSorghumjustMetabsSorghum)), 2)
  write.table(topMetabRounded, file = paste(fileForTopFeaturesjustMetabsSorghumSorghum,"sorghum_rounded", sep ="_"), row.names = FALSE)
  
  #output just the metabolites
  fileForTopFeaturesJustGenes = paste0(parentDir,"/drought_responsive_genes_only", paste(timeSample,collapse = '_'), ".txt")
  topFeatureInfoSorghumJustGenes = topFeatureInfoSorghum[topFeatureInfoSorghum %in% rownames(forPCATPMHeadersSorghum)]
  #write.table(topFeatureInfoSorghumJustGenes, file = fileForTopFeaturesJustGenes)
  
  responsiveInfoTableSorghum[rownames(responsiveInfoTableSorghum) == timeSample, colnames(responsiveInfoTableSorghum) == "responsive_metabolites"] = length(topFeatureInfoSorghumjustMetabsSorghum) 
  responsiveInfoTableSorghum[rownames(responsiveInfoTableSorghum) == timeSample, colnames(responsiveInfoTableSorghum) == "responsive_genes"] = length(topFeatureInfoSorghumJustGenes) 

  dev.off()
}




pdf("first_run_of_pipeline_complete_may_10th.pdf")
barplot(responsiveInfoTableSorghum$responsive_metabolites, names.arg = timePointsToParseSorghum, main = "Metabolites in Sorghum Impacted by Drought Stress Across Time Points as Determined by XCMS HILIC", ylab = "Number of Metabolites")

barplot(responsiveInfoTableSorghum$responsive_genes, names.arg = timePointsToParseSorghum, main = "Genes in Sorghum Impacted by Drought Stress Across Time Points as Determined by T-test", ylab = "Number of Genes")



#just HILIC data
barplot(plsDAInfoTablejustMetsSorghum$plsdaInfo, names.arg = timePointsToParseSorghum, main = "Error Rate of PLS-DA Model Each time Point with Sorghum HILIC Only", ylab = "Error Rate")

#just transcriptome
barplot(plsDAInfoTableJustRNA$plsdaInfo, names.arg = timePointsToParseSorghum, main = "Error Rate of best PLS-DA Model Each time Point with Sorghum RNA-Seq Only", ylab = "Error Rate")
dev.off()





barplot(plsDAInfoTablePermutations$plsdaInfo, names.arg = timePointsToParseSorghum, main = "Error Rate of best PLS-DA Model Each time Point with Sorghum HILIC And Transcriptome Permutation Test", ylab = "Error Rate")
plsDAInfoTablePermutations[rownames(plsDAInfoTablePermutations) == timeSample, colnames(plsDAInfoTablePermutations) == "bestCompNumber"] = which.min(perf.plsda$error.rate$overall)


#now, run everything as 


#have a table to store this info in

