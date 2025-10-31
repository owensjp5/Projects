# clear env
rm(list=ls())

# dataset
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150910

# load data
data <- read.table('/Users/jackowens/desktop/Projects/Transcriptomics/ref/GSE150910_gene-level_count_file.csv', header=TRUE, sep=',', row.names=1)

samples <- colnames(data)
genes <- rownames(data)

classes <- c()

# read class for each sample into classes
for(sample in samples){

    tmp <- unlist(strsplit(sample,'_'))
    classes <- append(classes, tmp[1])

}

# filter out chronic hypersensitivity pneumonitis(chp) data
filtered_data <- data[, classes %in% c('control','ipf')]
filtered_classes <- classes[classes %in% c('control','ipf')]

# save filtered data
saveRDS(filtered_data, file = "rds_objects/filtered_data.RDS")
saveRDS(filtered_classes, file = "rds_objects/filtered_classes.RDS")
