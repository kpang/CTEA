## Average the counts of each gene in all subtypes of each major cell type.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/expression_agg_4.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t")

    data_expression_1 <- rowMeans(data_expression[,2:ncol(data_expression)])

    data_expression_2 <- data.frame(data_expression[,1], data_expression_1)

    file_name <- paste("./", "/expression_agg_5.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    data_expression <- read.table(file_name, sep = "\t")

    data_expression_1 <- data_expression[data_expression[,2]>0,]

    file_name <- paste("./", "/high_expression.txt", sep = as.character(i))

    write.table(data_expression_1, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    file_name <- paste("./", "/gene_background.txt", sep = as.character(i))

    write.table(data_expression_1[,1], file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Rank all genes in each major cell type.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/high_expression.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t")

    data_rank <- rank(-data_expression[,2], ties.method = "average")

    data_rank <- 100*data_rank/length(data_rank)

    data_expression_1 <- data.frame(data_rank, data_expression)

    file_name <- paste("./", "/expression_rank.txt", sep = as.character(i))

    write.table(data_expression_1, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    data_expression <- read.table(file_name, sep = "\t")

    data_expression <- data_expression[ do.call(order, data_expression) ,]

    data_expression_1 <- data.frame(data_expression[,2:3], data_expression[,1])

    file_name <- paste("./", "/expression_rank_1.txt", sep = as.character(i))

    write.table(data_expression_1, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Extract the top 10% of genes with the highest expression levels in each major cell type.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/expression_rank_1.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t")

    data_expression_1 <- data_expression[data_expression[,3]<=10,]

    data_gene <- data_expression_1[,1]

    data_gene <- unique(data_gene)

    data_gene <- sort(data_gene)

    file_name <- paste("./", "/top_1000.txt", sep = as.character(i))

    write.table(data_gene, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Expression and percentage of the top 10% of genes.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/expression_rank_1.txt", sep = as.character(i))

    data_expression <- read.table(file_name, sep = "\t")

    data_expression_1 <- data_expression[data_expression[,3]<=10,]

    data_expression_2 <- data.frame(data_expression_1[,1], data_expression_1[,3], data_expression_1[,2])

    file_name <- paste("./", "/rank_top_1000.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Background genes for 30 major cell types.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

data_gene <- scan("./1/gene_background.txt", what = "character", sep = "\t")

for (i in 2:length(data_cell_type))
{
    file_name <- paste("./", "/gene_background.txt", sep = as.character(i))

    data_gene_1 <- scan(file_name, what = "character", sep = "\t")

    data_gene <- union(data_gene, data_gene_1)
}

data_gene <- sort(data_gene)

write.table(data_gene, file = "gene_background.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
