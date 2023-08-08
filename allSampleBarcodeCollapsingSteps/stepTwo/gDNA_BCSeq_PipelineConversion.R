#======================================================================================================================================
#edited by NJ/CJ to accommodate input from gDNA barcode pipeline
#======================================================================================================================================

#choose the analyzed folder which contains each of your individual samples after Step One
homeDirectory <- '/project/arjunrajlab/projectDirectory/experimentDirectory/analyzed/'
outputDirectory <- '/project/arjunrajlab/projectDirectory/experimentDirectory/Analysis/stepTwo/' #need to make directory before running

library(stringdist)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(spgs)

#import each sample table with barcodes and UMI counts into R for modificaiton
sampleFolders <- list.files(path = homeDirectory, pattern = '*UMICountsOnly', recursive = TRUE)
sampleTables <- list(rep("NA", length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(homeDirectory, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}
sampleNames <- dirname(dirname(sampleFolders))

#convert each sample table into the format expected by the side reaction pipeline
sampleSCTables <- list(rep("NA", length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleSCTable <-
    as_tibble(read.table(paste0(homeDirectory, sampleFolders[i]), stringsAsFactors = F)) %>% dplyr::rename(BC = V1, UMI = V2) %>%
    rowwise() %>% mutate(BC = reverseComplement(BC, case = "as is")) %>% ungroup() %>%
    mutate(BC50 = substring(BC, 1, 50), BC40 = substring(BC, 1, 40), BC30a = substring(BC, 1, 30), BC30b = substring(BC, 1, 30)) %>%
    mutate(sampleNum = sampleNames[i])
  sampleSCTables[[i]] <- sampleSCTable
}

sampleSCTableAll <- list()

#concatenate all sample tables into a single experiment table to feed into starcode
for(i in 1:length(sampleFolders)) {
  sampleSCTableTemp = sampleSCTables[[i]] %>% mutate(sampleNum = sampleNames[i])
  if (is.null(dim(sampleSCTableAll))) {sampleSCTableAll = sampleSCTableTemp}
  else {sampleSCTableAll = bind_rows(sampleSCTableAll, sampleSCTableTemp)}
}

#reformat the columns into the format expected by the side reaction pipeine
sampleSCTableAll <- sampleSCTableAll %>% add_column(cellID = "dummy", .before = "BC")
sampleSCTableAll <- sampleSCTableAll %>% relocate(UMI, .before = BC)
sampleSCTableAll$UMI <- as.character(sampleSCTableAll$UMI)

write.table(sampleSCTableAll, file= paste0(outputDirectory,'stepTwoCellIDUMIBarcodes.txt'),row.names=F,col.names=T,quote=F,sep="\t")
write.table(sampleSCTableAll[,4], file= paste0(outputDirectory,'stepTwoBarcodes50.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(sampleSCTableAll[,5], file= paste0(outputDirectory,'stepTwoBarcodes40.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(sampleSCTableAll[,6], file= paste0(outputDirectory,'stepTwoBarcodes30.txt'),row.names=F,col.names=F,quote=F,sep="\t")

#for checking LV distance histogram for an individual sample
barcodes = unique(sampleSCTableAll$BC)

set.seed(2059)
subsample1 = sample(barcodes, 1000)
subsample2 = sample(barcodes, 1000)
subsample3 = sample(barcodes, 1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

BarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic()

ggsave(filename = paste0(outputDirectory, "lvDistancePlot.pdf"), plot = BarcodesLvHistPlot, height = 5, width = 5)