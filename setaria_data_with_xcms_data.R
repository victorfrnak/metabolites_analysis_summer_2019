install.packages("VGAM")
library(VGAM)


fileSetaria = read.csv(file = "/Users/ahubbard/setaria_more_stringent_settings_minfrac_and_ppm_first_analysis.txt", sep = "\t", header = TRUE, comment.char = "",quote = "",stringsAsFactors = FALSE, na.strings = "NA")
namesForfileSetaria = paste(fileSetaria$mz, fileSetaria$rt, sep = "_")
rownames(fileSetaria) = namesForfileSetaria

colnames(fileSetaria) = gsub(".*_(.*)_.*", "\\1", colnames(fileSetaria))
colnames(fileSetaria) = gsub("\\.", "-", colnames(fileSetaria))


#change the names if needed from the xcms output
for(i in 1:length(colnames(fileSetaria)))
{
  
  
  colnames(fileSetaria)[i] = gsub(".y", "",as.character(colnames(fileSetaria)[i]))
  colnames(fileSetaria)[i] = gsub(".x", "",as.character(colnames(fileSetaria)[i]))
  #gsub('.x', colnames(fileSetaria)[i])
  toTest = strsplit(as.character(colnames(fileSetaria)[i]), "_")[[1]][length(strsplit(as.character(colnames(fileSetaria)[i]), "_")[[1]])]
  
  #print(i)
  #print(toTest)
  
  #if we need to, remove the last field
  if(grepl("T", as.character(toTest)) == FALSE)
  {
    toReplace = paste(strsplit(as.character(colnames(fileSetaria)[i]), "_")[[1]][-length(strsplit(as.character(colnames(fileSetaria)[i]), "_")[[1]])], collapse  = "_")
    print(toReplace)
    
    colnames(fileSetaria)[i] = toReplace
  }
  #remove the .x and .y
}


#vec of libs we'll want to pull out
#this is going to be metabolite data for which we
#also have transcriptome info

forPCATPMHeadersSetaria = read.table("/Users/ahubbard/drought_stress_joined_and_sorted_setaria_FPKM_tables_large_to_cluster.txt", sep = "\t", header = TRUE)
colnames(forPCATPMHeadersSetaria) = gsub(".*_P1_","",colnames(forPCATPMHeadersSetaria))
colnames(forPCATPMHeadersSetaria) = gsub(".*_P2_","",colnames(forPCATPMHeadersSetaria)) 

colnames(forPCATPMHeadersSetaria) = gsub("(.*)_.*", "\\1", colnames(forPCATPMHeadersSetaria))
colnames(forPCATPMHeadersSetaria) = gsub("(?<![0-9])0+", "", colnames(forPCATPMHeadersSetaria), perl = TRUE)
rownames(forPCATPMHeadersSetaria) = forPCATPMHeadersSetaria$Name
forPCATPMHeadersSetaria$Name = NULL



### read in the metadata ###
designSetupSetaria = read.table("/Users/ahubbard/Desktop/setaria_metabs_data_2018/hillic_data/metadata_complete_reformatted.txt", sep = "\t", header = TRUE)
designSetupSetariaNoOutliers = read.table("/Users/ahubbard/designSetupNoOutliers.txt", sep = "\t", header = TRUE)
designSetupSetariaNoOutliers$Tag = NULL
designSetupSetariaNoOutliers = designSetupSetariaNoOutliers
designSetupSetaria = designSetupSetaria[-1,]
designSetupSetaria = designSetupSetariaNoOutliers

#read in the setaria metadata
setariaMetadata = read.table("/Users/ahubbard/2018-10-30_Setaria_RNA-seq_sample_ID_coordinates_for_Allen_to_read_in.txt", sep = "\t", header = TRUE)
setariaMetadata$Just.Genotype = gsub(" ", "", setariaMetadata$Just.Genotype, fixed = TRUE) 
setariaMetadata$Just.Genotype = gsub(" ", "", setariaMetadata$Just.Genotype, fixed = TRUE) 
setariaMetadata$Just.treatment = gsub(" ", "", setariaMetadata$Just.treatment, fixed = TRUE) 


for (lib in 1:nrow(transcriptMetaSetaria))
{
  theLibInfo = transcriptMetaSetaria[rownames(transcriptMetaSetaria) == lib,]
  
  #paste together the time point, genotype, Treatment,rep
  #libGeno = theLibInfo[colnames(theLibInfo) == "Genotype.Treatment"] 
  libTreatment = theLibInfo[colnames(theLibInfo) == "Just.treatment"] 
  libDay = theLibInfo[colnames(theLibInfo) == "Day"]
  libAccession = theLibInfo[colnames(theLibInfo) == "Just.Genotype"]
  theReplicate = theLibInfo[colnames(theLibInfo) == "Rep"]
  theLibraryNumber = theLibInfo[colnames(theLibInfo) == "RNA.seq.sample.ID.on.tube"]
  
  
  #bind all of the info from the metadatas table together
  allTranscriptomeLibInfo = paste0(libDay, libAccession, libTreatment, theReplicate)
  allTranscriptomeLibInfo = str_replace_all(allTranscriptomeLibInfo, fixed(" "), "")
  vectorOfTranscriptInfo = c(vectorOfTranscriptInfo, allTranscriptomeLibInfo)
  RNAlibraryNumbers = c(RNAlibraryNumbers, theLibraryNumber)
  
  #replace with the correct format now
  
  #libGenoTreatment = gsub(" ","",libGenoTreatment)
  colnames(forPCATPMHeadersSetaria)[colnames(forPCATPMHeadersSetaria) == as.character(theLibraryNumber)] = allTranscriptomeLibInfo
}

timePointsToParse = as.vector(unique(designSetupSetaria$Time))

#don't cut out T0 this time
#timePointsToParse = timePointsToParse[-1]


#now, we'll join the transcriptome data
#read in the the transcriptome metadata
transcriptMetaSetaria = read.table(file= "/Users/ahubbard/2018-10-30_Setaria_RNA-seq_sample_ID_coordinates_for_Allen_to_read_in_names_edited.txt",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
forPCATPMHeadersSetaria = forPCATPMHeadersSetaria[rowMeans(forPCATPMHeadersSetaria) > 1,]

#convert RNA seq data to numeric 
forPCATPMHeadersSetaria[] <- lapply(forPCATPMHeadersSetaria, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

responsiveInfoTableSetaria <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParse)),y=0)) 
colnames(responsiveInfoTableSetaria) = c("responsive_metabolites","responsive_genes")
rownames(responsiveInfoTableSetaria) = timePointsToParse

 
vectorOfTranscriptInfo = vector()
RNAlibraryNumbers = vector()



forPCATPMHeadersSetaria = forPCATPMHeadersSetaria[rowMeans(forPCATPMHeadersSetaria) > 1,]




theCombinedGenesAndMetabNames = c(rownames(fileSetaria),rownames(forPCATPMHeadersSetaria))
infoOnEachMetaboliteSetaria <- data.frame(zeros(matrix(ncol = 2, nrow = length(theCombinedGenesAndMetabNames))))
colnames(infoOnEachMetaboliteSetaria) = c("Comparisons_drought_responsive","Which Time Points")
rownames(infoOnEachMetaboliteSetaria) = theCombinedGenesAndMetabNames



#parentDir = "/Users/ahubbard/dec5th_pipeline_totalII"
#parentDir = "/Users/ahubbard/april28th_pipeline_totalIII"
parentDir = getwd()
responsiveInfoTableSetaria <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParse)),y=0)) 
colnames(responsiveInfoTableSetaria) = c("responsive_metabolites","responsive_genes")
rownames(responsiveInfoTableSetaria) = timePointsToParse


plsDAInfoTableSetaria <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParse)),y=0)) 
colnames(plsDAInfoTableSetaria) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTableSetaria) = timePointsToParse


plsDAInfoTableSetariaJustMets <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParse)),y=0)) 
colnames(plsDAInfoTableSetariaJustMets) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTableSetariaJustMets) = timePointsToParse

plsDAInfoTableSetariaJustTrans <- data.frame(zeros(matrix(ncol = 2, nrow = length(timePointsToParse)),y=0)) 
colnames(plsDAInfoTableSetariaJustTrans) = c("plsdaInfo","bestCompNumber")
rownames(plsDAInfoTableSetariaJustTrans) = timePointsToParse


#put everything in here, this directory
setwd("/Users/ahubbard/xcms_metabolomics/setaria_may_13th_final/")
parentDir = getwd()

fileSetaria[fileSetaria == 0] <- NA
fileSetaria = na.omit(fileSetaria)

tableNoMissingSetaria = fileSetaria



#toPlot = fileSetaria[,colnames(fileSetaria) %in% intersect(colnames(fileSetaria), designSetupSetaria$Name)]
#metadatInfo = designSetupSetaria[designSetupSetaria$Name %in% intersect(colnames(fileSetaria), designSetupSetaria$Name),]


for( t in 1:(length(timePointsToParse)))
{
  
 
  #just make sure that we re-set the table each time
  tableNoMissingSetaria  = fileSetaria
  fileToJoinSetaria = t(tableNoMissingSetaria)
  timeSample = as.vector(unique(timePointsToParse)[t])
  
  pdf(paste(timeSample, "allPlotsMay13th.pdf", sep = "_"))
  
  
  theSubsetSetaria  = subset(designSetupSetaria, Time == timeSample)
  theDirForOutput = paste0(parentDir, "/",paste(timeSample,collapse = ""),"_full_report_all_one_file")
  mkdir(theDirForOutput)
  
  #because we'll take out the special characters later
  #but need to match with genotype now
  mkdir(paste0("/Users/ahubbard/PCA_testing_all_setaria_no_outliers/", timeSampleForDir))
  #file = theSubsetSetaria 
  
  
  #get a subset od all the libraries 
  #in this sample
  theSubsetSetaria  = subset(designSetupSetaria, Time %in% timeSample)
  
  
  #update the table of no missing values to be just the subset 
  tableNoMissingSetaria  = tableNoMissingSetaria[,colnames(tableNoMissingSetaria) %in% theSubsetSetaria$Name]
 
  
  tableNoMissingSetaria = as.data.frame(tableNoMissingSetaria)
  
  
  #seems that we have factors where we want numbers
  #but we don't want the NA's to be zeroes
  tableNoMissingSetaria[] <- lapply(tableNoMissingSetaria, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  
  forPCANoZeroesSetaria = tableNoMissingSetaria
  tableNoMissingSetaria = log(tableNoMissingSetaria + .0001)
  
  #rownames(forPCANoZeroesSetaria) = dataAndMetadataMerged$Name
  
  #name a table for the t-test 
  justMetsForTtest = tableNoMissingSetaria
  
  #before we do the tests, we want to convert to numeric 
  #run the pls-da and keep track
  #of the outputs
  
  XJustMetsSetaria = t(forPCANoZeroesSetaria)
  YJustMetsSetaria = rownames(t(forPCANoZeroesSetaria))
  
  YJustMetsSetaria = gsub(".*A10-1","",YJustMetsSetaria)
  YJustMetsSetaria = gsub(".*TB12-48","",YJustMetsSetaria)
  YJustMetsSetaria = gsub(".*TB12-201","",YJustMetsSetaria)

  YJustMetsSetaria = substr(YJustMetsSetaria,1,nchar(as.character(YJustMetsSetaria))-1)
  

  set.seed(2543) 
  plsdatestSetariaJustMets <- plsda(as.matrix(XJustMetsSetaria), as.factor(YJustMetsSetaria), ncomp = 2)
  plotIndiv(plsdatestSetariaJustMets, ind.names = TRUE, ellipse = TRUE, legend = TRUE, mainTitle = paste0("PLS-DA at Combined Metabolome + Transcriptome", timePointsToParse))
  
  
  set.seed(2543) 
  perf.plsdaSetariaJustMets  <- perf(plsdatestSetariaJustMets, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = 10)
  #plot(perf.plsdaSetariaJustMets, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal", mainTitle = paste0("PLS-DA at Combined Metabolome + Transcriptome", timePointsToParse))
  
  myPLSDAErrorSetariaJustMets = perf.plsdaSetariaJustMets$error.rate$overall[rownames(perf.plsdaSetariaJustMets$error.rate$overall) == "comp 2",colnames(perf.plsdaSetariaJustMets$error.rate$overall) == "centroids.dist"]
  
  plsDAInfoTableSetariaJustMets[rownames(plsDAInfoTableSetariaJustMets) == timeSample, colnames(plsDAInfoTableSetariaJustMets) == "plsdaInfo"] = myPLSDAErrorSetariaJustMets
  
  
  #do the combined transcriptome and metabolome, now
  
  
  
  treatmentPCA = paste0(theDirForOutput, "/allTimepointByTreatment.pdf")

  toMergeWithMetabsSetaria = forPCATPMHeadersSetaria[,colnames(forPCATPMHeadersSetaria) %in% colnames(tableNoMissingSetaria)]
  
  
  
  #do just the transcriptome now
  XjustTransSetaria = t(toMergeWithMetabsSetaria)
  YjustTransSetaria = rownames(t(toMergeWithMetabsSetaria))
  
  YjustTransSetaria = gsub(".*A10-1","",YjustTransSetaria)
  YjustTransSetaria = gsub(".*TB12-48","",YjustTransSetaria)
  YjustTransSetaria = gsub(".*TB12-201","",YjustTransSetaria)
  #YjustTransSetaria = gsub("Rep_.*","",Y)
  
  
  YjustTransSetaria = substr(YjustTransSetaria,1,nchar(as.character(YjustTransSetaria))-1)
  
  
  set.seed(2543) 
  plsdatestSetariajustTrans <- plsda(as.matrix(XjustTransSetaria), as.factor(YjustTransSetaria), ncomp = 2)
  plotIndiv(plsdatestSetariajustTrans, ind.names = TRUE, ellipse = TRUE, legend = TRUE, mainTitle = paste0("PLS-DA at Combined Metabolome + Transcriptome", timePointsToParse))
  
  
  set.seed(2543) 
  perf.plsdaSetariajustTrans  <- perf(plsdatestSetariajustTrans, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = 10)
  plot(perf.plsdaSetariajustTrans , col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal", mainTitle = paste0("PLS-DA at Combined Metabolome + Transcriptome", timePointsToParse))
  
  myPLSDAErrorSetariaJustTrans = perf.plsdaSetariajustTrans$error.rate$overall[rownames(perf.plsdaSetariajustTrans$error.rate$overall) == "comp 2",colnames(perf.plsdaSetariajustTrans$error.rate$overall) == "centroids.dist"]
  
  
  plsDAInfoTableSetariaJustTrans[rownames(plsDAInfoTableSetariaJustTrans) == timeSample, colnames(plsDAInfoTableSetariaJustTrans) == "plsdaInfo"] = myPLSDAErrorSetariaJustTrans
  
  
  
  #focus on only the compounds that are drought responsive
  #at one or more time points
  topFeatureInfo = vector()  
  
  #focus only on the columns that are both in the transcriptome and metabolome
  tableNoMissingSetaria = tableNoMissingSetaria[,colnames(tableNoMissingSetaria) %in% colnames(forPCATPMHeadersSetaria)]
  #toMergeWithMetabsSetaria = forPCATPMHeadersSetaria[,colnames(forPCATPMHeadersSetaria) %in% colnames(tableNoMissingSetaria)]
  
  justMetsForTtest = tableNoMissingSetaria
  #tableNoMissingSetaria = rbind(tableNoMissingSetaria, toMergeWithMetabsSetaria)
  
  #in subset in order to do t-test
  #tableNoMissingSetaria = justMetsForTtest
  
  
  #pooling the genotypes together into just control an drought stress
  
    for (a in 1:nrow(tableNoMissingSetaria))
    {
  
     theCompound = rownames(tableNoMissingSetaria)[a]
  
      timePoint = timeSample 
  
      vec40Percent = vector()
      vec100Percent = vector()
  
      theSubsetSetaria  = subset(designSetupSetaria, Time %in% timeSample)
      for (k in 1: (length(unique(designSetupSetaria$Treatment))))
      {
  
        treatmentSample = unique(designSetupSetaria$Treatment)[k]  
  
        theSubsetSetariaII = theSubsetSetaria[theSubset$Treatment == treatmentSample,]
        theSamplesToImputeSetaria = theSubsetSetariaII$Name
  
  
         ### get the info for the metabolite ###
        forTtest = tableNoMissingSetaria[rownames(tableNoMissingSetaria) == theCompound,intersect(colnames(tableNoMissingSetaria), theSamplesToImputeSetaria)]
  
  
        if(treatmentSample == 100)
        {
          vec100Percent = forTtest
          #print("hello")
        }
  
        if(treatmentSample == 40)
        {
          vec40Percent = forTtest
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
        infoOnEachMetaboliteSetaria[rownames(infoOnEachMetaboliteSetaria) == theCompound, 1] =  infoOnEachMetaboliteSetaria[rownames(infoOnEachMetaboliteSetaria) == theCompound, 1] + 1
        tTestResults = theTest$estimate
  
        if(tTestResults[1] > tTestResults[2])
        {
          resultsForTable = paste0("upInDrought",timeSample) 
        }
        if(tTestResults[2] > tTestResults[1])
       {
         resultsForTable = paste0("downInDrought",timeSample) 
        }
  
        infoOnEachMetaboliteSetaria[rownames(infoOnEachMetaboliteSetaria) == theCompound, 2] =  paste(infoOnEachMetaboliteSetaria[rownames(infoOnEachMetaboliteSetaria) == theCompound, 2], resultsForTable, collapse = ",")
        topFeatureInfo = c(topFeatureInfo, theCompound)
  
      #we'll also want to do by genotype, I assume 
      }
    }
  
  
  #output just the genes
  fileForTopFeaturesJustMetabs = paste0(parentDir,"/drought_responsive_metabolites_only_setaria_may12th", paste(timeSample,collapse = ''), ".txt")
  topFeatureInfoJustMetabs = topFeatureInfo[topFeatureInfo %in% rownames(fileSetaria)]
  
  #split into the numbers and round up
  topFeatureInfoJustMetabs = gsub("(.*)_.*", "\\1", topFeatureInfoJustMetabs)
  
  #write.table(topFeatureInfoJustMetabs, file = fileForTopFeaturesJustMetabs)
  topMetabRounded = round(as.numeric(gsub("(.*)_.*", "\\1", topFeatureInfoJustMetabs)), 2)
  write.table(topMetabRounded, file = paste(fileForTopFeaturesJustMetabs,"setaria_rounded", sep ="_"), row.names = FALSE)
  
  #output just the metabolites
  fileForTopFeaturesJustGenes = paste0(parentDir,"/drought_responsive_genes_only", paste(timeSample,collapse = ''), ".txt")
  topFeatureInfoJustGenes = topFeatureInfo[topFeatureInfo %in% rownames(forPCATPMHeadersSetaria)]
  #write.table(topFeatureInfoJustGenes, file = fileForTopFeaturesJustGenes)
  
  responsiveInfoTableSetaria[rownames(responsiveInfoTableSetaria) == timeSample, colnames(responsiveInfoTableSetaria) == "responsive_metabolites"] = length(topFeatureInfoJustMetabs) 
  responsiveInfoTableSetaria[rownames(responsiveInfoTableSetaria) == timeSample, colnames(responsiveInfoTableSetaria) == "responsive_genes"] = length(topFeatureInfoJustGenes) 
  
  dev.off()
}

barplot(plsDAInfoTableSetaria$plsdaInfo, names.arg = timePointsToParse, main = "Error Rate of best PLS-DA Model Each time Point with Setaria Combined", ylab = "Error Rate")

barplot(plsDAInfoTableSetaria$plsdaInfo, names.arg = timePointsToParse, main = "Error Rate of best PLS-DA Model Each time Point with Setaria Combined", ylab = "Error Rate")



barplot(plsDAInfoTableSetariaJustMets$plsdaInfo, names.arg = timePointsToParse, main = "Error Rate of PLS-DA Model Each time Point with Setaria HILIC Only", ylab = "Number of Metabolites")



barplot(responsiveInfoTableSetaria$responsive_metabolites, names.arg = timePointsToParse, main = "Metabolites in Setaria Impacted by Drought Stress Across Time Points as Determined by XCMS HILIC", ylab = "Number of Metabolites")




barplot(responsiveInfoTableSetaria$responsive_genes, names.arg = timePointsToParse, main = "Metabolites in Setaria Impacted by Drought Stress Across Time Points as Determined by XCMS HILIC", ylab = "Number of Metabolites")
