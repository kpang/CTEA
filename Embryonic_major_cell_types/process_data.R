## Read the aggregate cell-cluster-level scRNA-seq data of 794 cell clusters from 165 cell subtypes and 22 major cell types in the embryonic mouse brain.

library(loomR)

## The file "dev_all.agg.loom" contains the aggregate scRNA-seq data of 22 major cell types in the embryonic mouse brain. This data file can be downloaded from http://mousebrain.org/development/downloads.

lfile <- connect(filename = "dev_all.agg.loom", mode = "r", skip.validate = TRUE)

data_expression <- lfile[["matrix"]][, ]

data_expression <- t(data_expression)

data_gene <- lfile[["row_attrs/Gene"]][]

data_expression_1 <- data.frame(data_gene, data_expression)

write.table(data_expression_1, file = "expression_agg.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(data_gene, file = "gene_agg.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_sample_1 <- lfile[["col_attrs/Clusters"]][]

data_sample_2 <- lfile[["col_attrs/ClusterName"]][]

data_sample_3 <- lfile[["col_attrs/Subclass"]][]

data_sample_4 <- lfile[["col_attrs/Class"]][]

data_sample <- data.frame(data_sample_1, data_sample_2)

write.table(data_sample, file = "sample_id_agg.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_sample <- data.frame(data_sample_1, data_sample_3)

write.table(data_sample, file = "sample_id_agg_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_sample <- data.frame(data_sample_1, data_sample_4)

write.table(data_sample, file = "sample_id_agg_2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

lfile$close_all()

## Normalize the scRNA-seq data through dividing the count of each gene in a cell cluster by the total count of the cell cluster and multiplying it by 100,000.

data_expression <- read.table("expression_agg.txt", sep = "\t")

data_col <- colSums(data_expression[,2:ncol(data_expression)])

for (i in 1:(ncol(data_expression)-1))
{
    data_expression[,i+1] <- 100000*data_expression[,i+1]/data_col[i]
}

write.table(data_expression, file = "expression_agg_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Average expression for duplicated mouse genes.

data_expression <- read.table("expression_agg_1.txt", sep = "\t")

data_gene <- data_expression[,1]

data_gene_1 <- unique(data_gene)

data_gene_2 <- data_gene[duplicated(data_gene)]

data_gene_2 <- unique(data_gene_2)

data_gene_2 <- sort(data_gene_2)

data_gene_3 <- setdiff(data_gene_1, data_gene_2)

data_gene_3 <- sort(data_gene_3)

data_expression_1 <- merge(data_expression, data_gene_2, by = 1)

data_expression_2 <- merge(data_expression, data_gene_3, by = 1)

write.table(data_expression_2, file = "expression_agg_2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

for (i in 1:length(data_gene_2))
{
    data_expression_3 <- data_expression_1[data_expression_1[,1]==data_gene_2[i],]

    data_expression_4 <- colMeans(data_expression_3[,2:ncol(data_expression_3)])

    data_expression_5 <- c(as.character(data_gene_2[i]), data_expression_4)

    write.table(t(data_expression_5), file = "expression_agg_2.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Remove genes with count of zero across all cell clusters.

data_expression <- read.table("expression_agg_2.txt", sep = "\t")

data_row <- rowSums(data_expression[,2:ncol(data_expression)]>0)

data_expression_1 <- data_expression[data_row>=1,]

write.table(data_expression_1, file = "expression_agg_3.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_gene <- data_expression_1[,1]

write.table(data_gene, file = "gene_agg_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Annotation data of subtypes.

data_sample <- read.table("sample_id_agg_2.txt", sep = "\t", colClasses = "character")

data_sample_1 <- 1:nrow(data_sample)

data_sample_2 <- data.frame(data_sample[,2], data_sample_1)

data_sample_2[379,1] <- "Neuroblast"

write.table(data_sample_2, file = "sample_id_agg_3.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Names of 22 major cell types.

data_sample <- read.table("sample_id_agg_3.txt", sep = "\t")

data_cell_type <- data_sample[,1]

data_cell_type <- unique(data_cell_type)

data_cell_type <- sort(data_cell_type)

write.table(data_cell_type, file = "sample_id_agg_4.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Divide the expression data into 22 major cell types.

data_sample <- read.table("sample_id_agg_3.txt", sep = "\t", colClasses = "character")

## The file "sample_id_agg_5.txt" is derived from the file "sample_id_agg_4.txt" with one major cell type removed.

data_cell_type <- scan("sample_id_agg_5.txt", what = "character", sep = "\t")

data_expression <- read.table("expression_agg_3.txt", sep = "\t")

data_gene <- data_expression[,1]

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_position <- data_position+1

    data_expression_1 <- data_expression[,data_position]

    data_expression_2 <- data.frame(data_gene, data_expression_1)

    dir.create(as.character(i))

    file_name <- paste("./", "/expression_agg_2.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Divide the annotation data into 22 major cell types.

data_sample <- read.table("sample_id_agg_3.txt", sep = "\t", colClasses = "character")

data_cell_type <- scan("sample_id_agg_5.txt", what = "character", sep = "\t")

data_sample_1 <- read.table("sample_id_agg.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_sample_2 <- data_sample_1[data_position,]

    file_name <- paste("./", "/sample_id_agg.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data_sample_1 <- read.table("sample_id_agg_1.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_sample_2 <- data_sample_1[data_position,]

    file_name <- paste("./", "/sample_id_agg_1.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data_sample_1 <- read.table("sample_id_agg_2.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_sample_2 <- data_sample_1[data_position,]

    file_name <- paste("./", "/sample_id_agg_2.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Map mouse gene symbol to human gene entrezid.

data_cell_type <- scan("sample_id_agg_5.txt", what = "character", sep = "\t")

## The file "mouse_symbol2human_entrezid.txt" contains gene ID mapping between mouse gene symbol and human gene entrezid. This file was downloaded from the Ensembl database (http://www.ensembl.org).

data_gene <- read.table("mouse_symbol2human_entrezid.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/expression_agg_2.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t", colClasses = "character")

    data_expression_1 <- merge(data_gene, data_expression, by = 1)

    data_expression_2 <- data_expression_1[,2:ncol(data_expression_1)]

    data_expression_2 <- data_expression_2[ do.call(order, data_expression_2) ,]

    file_name <- paste("./", "/expression_agg_3.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Average expression for duplicated human genes.

data_cell_type <- scan("sample_id_agg_5.txt", what = "character", sep = "\t")

for (i in 1:20)
{
    file_name <- paste("./", "/expression_agg_3.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t")

    data_gene <- data_expression[,1]

    data_gene_1 <- unique(data_gene)

    data_gene_2 <- data_gene[duplicated(data_gene)]

    data_gene_2 <- unique(data_gene_2)

    data_gene_2 <- sort(data_gene_2)

    data_gene_3 <- setdiff(data_gene_1, data_gene_2)

    data_gene_3 <- sort(data_gene_3)

    data_expression_1 <- merge(data_expression, data_gene_2, by = 1)

    data_expression_2 <- merge(data_expression, data_gene_3, by = 1)

    file_name <- paste("./", "/expression_agg_4.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    for (j in 1:length(data_gene_2))
    {
        data_expression_3 <- data_expression_1[data_expression_1[,1]==data_gene_2[j],]

        data_expression_4 <- colMeans(data_expression_3[,2:ncol(data_expression_3)])

        data_expression_5 <- c(as.character(data_gene_2[j]), data_expression_4)

        write.table(t(data_expression_5), file = file_name, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
}

i <- 21

file_name <- paste("./", "/expression_agg_3.txt", sep = as.character(i))

data_expression <- read.table(file_name, sep = "\t")

data_gene <- data_expression[,1]

data_gene_1 <- unique(data_gene)

data_gene_2 <- data_gene[duplicated(data_gene)]

data_gene_2 <- unique(data_gene_2)

data_gene_2 <- sort(data_gene_2)

data_gene_3 <- setdiff(data_gene_1, data_gene_2)

data_gene_3 <- sort(data_gene_3)

data_expression_1 <- merge(data_expression, data_gene_2, by = 1)

data_expression_2 <- merge(data_expression, data_gene_3, by = 1)

file_name <- paste("./", "/expression_agg_4.txt", sep = as.character(i))

write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

for (j in 1:length(data_gene_2))
{
    data_expression_3 <- data_expression_1[data_expression_1[,1]==data_gene_2[j],]

    data_expression_4 <- mean(data_expression_3[,2:ncol(data_expression_3)])

    data_expression_5 <- c(as.character(data_gene_2[j]), data_expression_4)

    write.table(t(data_expression_5), file = file_name, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

i <- 22

file_name <- paste("./", "/expression_agg_3.txt", sep = as.character(i))

data_expression <- read.table(file_name, sep = "\t")

data_gene <- data_expression[,1]

data_gene_1 <- unique(data_gene)

data_gene_2 <- data_gene[duplicated(data_gene)]

data_gene_2 <- unique(data_gene_2)

data_gene_2 <- sort(data_gene_2)

data_gene_3 <- setdiff(data_gene_1, data_gene_2)

data_gene_3 <- sort(data_gene_3)

data_expression_1 <- merge(data_expression, data_gene_2, by = 1)

data_expression_2 <- merge(data_expression, data_gene_3, by = 1)

file_name <- paste("./", "/expression_agg_4.txt", sep = as.character(i))

write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

for (j in 1:length(data_gene_2))
{
    data_expression_3 <- data_expression_1[data_expression_1[,1]==data_gene_2[j],]

    data_expression_4 <- colMeans(data_expression_3[,2:ncol(data_expression_3)])

    data_expression_5 <- c(as.character(data_gene_2[j]), data_expression_4)

    write.table(t(data_expression_5), file = file_name, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
