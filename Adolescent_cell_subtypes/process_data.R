## Read the expression data of subtypes of adolescent telencephalon projecting excitatory neurons.

## The file "sample_id_agg_3.txt" contains names of subtypes of adolescent telencephalon projecting excitatory neurons.

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

data_expression <- read.table("../Adolescent_major_cell_types/26/expression_agg_3.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    dir.create(as.character(i))

    data_expression_1 <- data.frame(data_expression[,1], data_expression[,i+1])

    file_name <- paste("./", "/expression_agg_4.txt", sep = as.character(i))

    write.table(data_expression_1, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Read the annotation data of subtypes of adolescent telencephalon projecting excitatory neurons.

data_sample <- read.table("../Adolescent_major_cell_types/26/sample_id_agg.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/sample_id_agg.txt", sep = as.character(i))

    write.table(data_sample[i,], file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data_sample <- read.table("../Adolescent_major_cell_types/26/sample_id_agg_3.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    file_name <- paste("./", "/sample_id_agg_3.txt", sep = as.character(i))

    write.table(data_sample[i,], file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data1 <- read.table("../Adolescent_major_cell_types/26/sample_id_agg.txt", sep = "\t", colClasses = "character")

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

data2 <- data.frame(data1, data_cell_type)

write.table(data2, file = "sample_id_agg_4.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
