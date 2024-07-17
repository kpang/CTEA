## 30 sets of ASD and DD genes expressed in any of 22 major cell types of embryonic mouse brain.

data_background <- scan("gene_background.txt", what = "character", sep = "\t")

## The file "disease_set_1.txt" contains names of 30 ASD and DD gene set files. The genes in these gene set files can be obtained from Supplementary Table 1.

data_disease_set <- scan("disease_set_1.txt", what = "character", sep = "\t")

dir.create("disease_set")

for (i in 1:length(data_disease_set))
{
    file_name <- paste(data_disease_set[i], ".txt", sep = "")

    data_gene <- scan(file_name, what = "character", sep = "\t")

    data_gene_1 <- intersect(data_background, data_gene)

    data_gene_1 <- sort(data_gene_1)

    file_name_1 <- paste("./disease_set/", ".txt", sep = data_disease_set[i])

    write.table(data_gene_1, file = file_name_1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Background genes for each major cell type.

## The file "cell_type_1.txt" contains orders in which the 22 major cell types are shown in the plot.

data_cell_type <- scan("cell_type_1.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/gene_background.txt", sep = data_cell_type[i])

    data_gene <- scan(file_name, what = "character", sep = "\t")

    file_name_1 <- paste("background_gene_", ".txt", sep = as.character(i))

    write.table(data_gene, file = file_name_1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Top 10% of genes with the highest expression levels in each major cell type.

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/top_1000.txt", sep = data_cell_type[i])

    data_gene <- scan(file_name, what = "character", sep = "\t")

    file_name_1 <- paste("function_gene_", ".txt", sep = as.character(i))

    write.table(data_gene, file = file_name_1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## 30 sets of ASD and DD genes expressed in any major cell types.

data_disease_set <- scan("disease_set_1.txt", what = "character", sep = "\t")

for (i in 1:length(data_disease_set))
{
    file_name <- paste("./disease_set/", ".txt", sep = data_disease_set[i])

    data_gene <- scan(file_name, what = "character", sep = "\t")

    file_name_1 <- paste("disease_gene_", ".txt", sep = as.character(i))

    write.table(data_gene, file = file_name_1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Size of 30 sets of ASD and DD genes expressed in any major cell types.

for (i in 1:length(data_disease_set))
{
    file_name <- paste("./disease_set/", ".txt", sep = data_disease_set[i])

    data_gene <- scan(file_name, what = "character", sep = "\t")

    data_size <- length(data_gene)

    write.table(data_size, file = "set_size_1.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Calculate overlap between 30 ASD and DD gene sets and the top 10% of genes highly expressed in each of 30 major cell types.

data_cell_type <- scan("cell_type_1.txt", what = "character", sep = "\t")

data_disease_set <- scan("disease_set_1.txt", what = "character", sep = "\t")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("background_gene_", ".txt", sep = as.character(i))

    data_background <- scan(file_name, what = "character", sep = "\t")

    file_name <- paste("function_gene_", ".txt", sep = as.character(i))

    data_function_gene <- scan(file_name, what = "character", sep = "\t")

    file_name_1 <- paste("overlap_gene_", ".txt", sep = as.character(i))

    for (j in 1:length(data_disease_set))
    {
        file_name <- paste("disease_gene_", ".txt", sep = as.character(j))

        data_disease_gene <- scan(file_name, what = "character", sep = "\t")

        data_disease_gene <- intersect(data_disease_gene, data_background)

        data_intersect <- intersect(data_disease_gene, data_function_gene)

        data_union <- union(data_disease_gene, data_function_gene)

        data_intersect_union <- matrix(data = NA, nrow = 2, ncol = 2)

        data_intersect_union[1,1] <- length(data_intersect)

        data_intersect_union[1,2] <- length(data_disease_gene)-length(data_intersect)

        data_intersect_union[2,1] <- length(data_function_gene)-length(data_intersect)

        data_intersect_union[2,2] <- length(data_background)-length(data_union)

        data_test <- fisher.test(data_intersect_union, alternative = "greater")

        data_estimate <- data_test$estimate

        data_p <- data_test$p.value

        data_n <- length(data_intersect)

        data_n_1 <- length(data_function_gene)

        data_ratio_1 <- length(data_intersect)/length(data_function_gene)

        data_n_2 <- length(data_disease_gene)

        data_ratio_2 <- length(data_intersect)/length(data_disease_gene)

        data_overlap <- c(data_estimate, data_p, data_n, data_ratio_1, data_n_1, data_ratio_2, data_n_2)

        write.table(t(data_overlap), file = file_name_1, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
}

## Odds ratio and P-value of the overlap.

data_cell_type <- scan("cell_type_1.txt", what = "character", sep = "\t")

i <- 1

file_name <- paste("overlap_gene_", ".txt", sep = as.character(i))

data_overlap <- read.table(file_name, sep = "\t", colClasses = "character")

data_overlap_or <- data_overlap[,1]

data_overlap_p <- data_overlap[,2]

for (i in 2:length(data_cell_type))
{
    file_name <- paste("overlap_gene_", ".txt", sep = as.character(i))

    data_overlap <- read.table(file_name, sep = "\t", colClasses = "character")

    data_overlap_or <- data.frame(data_overlap_or, data_overlap[,1])

    data_overlap_p <- data.frame(data_overlap_p, data_overlap[,2])
}

write.table(data_overlap_or, file = "overlap_or.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(data_overlap_p, file = "overlap_p.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_overlap_or <- read.table("overlap_or.txt", sep = "\t")

data_overlap_or <- t(data_overlap_or)

write.table(data_overlap_or, file = "overlap_or_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_overlap_or <- read.table("overlap_or_1.txt", sep = "\t")

data_overlap_or <- signif(data_overlap_or, digits = 3)

write.table(data_overlap_or, file = "overlap_or_2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_overlap_p <- read.table("overlap_p.txt", sep = "\t")

data_overlap_p <- -log10(data_overlap_p)

data_overlap_p <- t(data_overlap_p)

write.table(data_overlap_p, file = "overlap_p_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_overlap_p <- read.table("overlap_p_1.txt", sep = "\t")

data_overlap_p <- signif(data_overlap_p, digits = 3)

write.table(data_overlap_p, file = "overlap_p_2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
