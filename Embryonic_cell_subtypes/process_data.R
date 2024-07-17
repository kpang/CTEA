## Read the annotation data of subtypes of embryonic neurons.

data_sample <- read.table("../Embryonic_major_cell_types/15/sample_id_agg_1.txt", sep = "\t", colClasses = "character")

data_position <- 1:nrow(data_sample)

data_sample_1 <- data.frame(data_sample[,2], data_position)

write.table(data_sample_1, file = "sample_id_agg_4.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Names of embryonic neuron subtypes.

data_sample <- read.table("sample_id_agg_4.txt", sep = "\t", colClasses = "character")

data_cell_type <- data_sample[,1]

data_cell_type <- unique(data_cell_type)

data_cell_type <- sort(data_cell_type)

write.table(data_cell_type, file = "sample_id_agg_3.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Number of cell clusters in each embryonic neuron subtype.

data_sample <- read.table("sample_id_agg_4.txt", sep = "\t", colClasses = "character")

data_sample <- data_sample[,1]

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

data_n <- numeric(length = length(data_cell_type))

for (i in 1:length(data_cell_type))
{
    data_sample_1 <- data_sample[data_sample==data_cell_type[i]]

    data_n[i] <- length(data_sample_1)
}

data_cell_type_1 <- data.frame(data_cell_type, data_n)

write.table(data_cell_type_1, file = "sample_id_agg_5.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Divide the expression data of embryonic neurons into cell subtypes.

data_sample <- read.table("sample_id_agg_4.txt", sep = "\t", colClasses = "character")

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

data_expression <- read.table("../Embryonic_major_cell_types/15/expression_agg_4.txt", sep = "\t")

data_gene <- data_expression[,1]

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_position <- data_position+1

    data_expression_1 <- data_expression[,data_position]

    data_expression_2 <- data.frame(data_gene, data_expression_1)

    dir.create(as.character(i))

    file_name <- paste("./", "/expression_agg_3.txt", sep = as.character(i))

    write.table(data_expression_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

## Divide the annotation data of embryonic neurons into cell subtypes.

data_sample <- read.table("sample_id_agg_4.txt", sep = "\t", colClasses = "character")

data_cell_type <- scan("sample_id_agg_3.txt", what = "character", sep = "\t")

data_sample_1 <- read.table("../Embryonic_major_cell_types/15/sample_id_agg.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_sample_2 <- data_sample_1[data_position,]

    file_name <- paste("./", "/sample_id_agg.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data_sample_1 <- read.table("../Embryonic_major_cell_types/15/sample_id_agg_1.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_sample_2 <- data_sample_1[data_position,]

    file_name <- paste("./", "/sample_id_agg_1.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data_sample_1 <- read.table("../Embryonic_major_cell_types/15/sample_id_agg_2.txt", sep = "\t", colClasses = "character")

for (i in 1:length(data_cell_type))
{
    data_position <- data_sample[data_sample[,1]==data_cell_type[i],2]

    data_position <- as.numeric(data_position)

    data_sample_2 <- data_sample_1[data_position,]

    file_name <- paste("./", "/sample_id_agg_2.txt", sep = as.character(i))

    write.table(data_sample_2, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
