## Read the aggregate cell-type-level scRNA-seq data of 253 cell subtypes from 30 major cell types in the adolescent mouse nervous system.

library(loomR)

## The file "aggregate_data_file_name.txt" contains names of 30 aggregate scRNA-seq data files of 30 major cell types in the adolescent mouse nervous system. These data files can be downloaded from http://mousebrain.org/adolescent/loomfiles_level_L6.

data_file_name <- scan("aggregate_data_file_name.txt", what = "character", sep = "\t")

for (i in 1:length(data_file_name))
{
    dir.create(as.character(i))

    lfile <- connect(filename = data_file_name[i], mode = "r")

    data_expression <- lfile[["matrix"]][, ]

    data_expression <- t(data_expression)

    data_gene <- lfile[["row_attrs/Gene"]][]

    data_expression_1 <- data.frame(data_gene, data_expression)

    file_name <- paste("./", "/expression_agg.txt", sep = as.character(i))

    write.table(data_expression_1, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    file_name <- paste("./", "/gene_agg.txt", sep = as.character(i))

    write.table(data_gene, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    data_sample_1 <- lfile[["col_attrs/Clusters"]][]

    data_sample_2 <- lfile[["col_attrs/Probable_location"]][]

    data_sample <- data.frame(data_sample_1, data_sample_2)

    file_name <- paste("./", "/sample_id_agg.txt", sep = as.character(i))

    write.table(data_sample, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    file_name <- paste("./", "/sample_id_agg_3.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    lfile$close_all()
}

## Normalize the scRNA-seq data through dividing the count of each gene in a cell subtype by the total count of the cell subtype and multiplying it by 100,000.

## The file "sample_id_agg_3.txt" contains names of 30 major cell types of adolescent mouse nervous system.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/expression_agg.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t")

    data_col <- colSums(data_expression[,2:ncol(data_expression)])

    for (j in 1:(ncol(data_expression)-1))
    {
        data_expression[,j+1] <- 100000*data_expression[,j+1]/data_col[j]
    }

    file_name <- paste("./", "/expression_agg_1.txt", sep = as.character(i))

    write.table(data_expression, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Average expression for duplicated mouse genes.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/expression_agg_1.txt", sep = as.character(i))

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

    file_name <- paste("./", "/expression_agg_2.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    for (j in 1:length(data_gene_2))
    {
        data_expression_3 <- data_expression_1[data_expression_1[,1]==data_gene_2[j],]

        data_expression_4 <- colMeans(data_expression_3[,2:ncol(data_expression_3)])

        data_expression_5 <- c(as.character(data_gene_2[j]), data_expression_4)

        write.table(t(data_expression_5), file = file_name, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
}

## Map mouse gene symbol to human gene entrezid.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

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

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
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
