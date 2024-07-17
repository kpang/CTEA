## Prepare data for overlap plot.

library(gplots)

data_cell_type <- scan("cell_type_1.txt", sep = "\t")

data_overlap_p <- read.table("overlap_p_2.txt", sep = "\t")

data_overlap_or <- read.table("overlap_or_2.txt", sep = "\t")

for (i in 1:ncol(data_overlap_p))
{
    data_cell_type_1 <- 1:length(data_cell_type)

    data_cell_type_2 <- rep(ncol(data_overlap_p)-i+1, length(data_cell_type))

    data_overlap_p_1 <- data_overlap_p[,i]

    data_overlap_or_1 <- data_overlap_or[,i]

    data_score <- data.frame(data_cell_type_1, data_cell_type_2, data_overlap_p_1, data_overlap_or_1, data_overlap_or_1)

    write.table(data_score, file = "cell_subtype.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data_score <- read.table("cell_subtype.txt", sep = "\t")

data_color <- colorpanel(6, "white", "red")

data_color <- gsub("#", "", data_color)

data_score[data_score[,4]<1,5] <- data_color[1]

data_score[(data_score[,4]>=1)&(data_score[,4]<2),5] <- data_color[2]

data_score[(data_score[,4]>=2)&(data_score[,4]<3),5] <- data_color[3]

data_score[(data_score[,4]>=3)&(data_score[,4]<4),5] <- data_color[4]

data_score[(data_score[,4]>=4)&(data_score[,4]<5),5] <- data_color[5]

data_score[data_score[,4]>=5,5] <- data_color[6]

write.table(data_score, file = "cell_subtype_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Plot overlap between 30 ASD and DD gene sets and the top 10% of genes highly expressed in each subtype of adolescent telencephalon projecting excitatory neurons.

data_cell_type <- scan("cell_type_1.txt", sep = "\t")

## The file "gene_set_1.txt" is the file "sample_id_agg_3.txt" with small changes.

data_function_set <- scan("gene_set_1.txt", what = "character", sep = "\t")

data_function_set <- data_function_set[data_cell_type]

data_disease_set <- scan("disease_set_2.txt", what = "character", sep = "\t")

data_set_size <- scan("set_size_1.txt", what = "character", sep = "\t")

data_disease_set <- paste(data_disease_set, data_set_size, sep = " (")

data_disease_set <- paste(data_disease_set, ")", sep = "")

data_width <- length(data_cell_type)*0.64

data_score <- read.table("cell_subtype_1.txt", sep = "\t")

data_color <- paste("#", as.character(data_score[,5]), sep = "")

pdf(file = "overlap_1.pdf", width = data_width, height = 15.5)

par(mar = c(21.5, 30, 2, 0))

par(mgp = c(3, 1, 0))

par(cex.axis = 1.5)

par(cex.lab = 1.5)

par(cex.main = 2)

par(bg = 0)

plot(data_score[,1], data_score[,2], type = "p", pch = 19, xlim = c(1, max(data_score[,1])), ylim = c(0, max(data_score[,2])), main = "", xlab = "", ylab = "", col = data_color, cex = (data_score[,3]+5)/5, axes = FALSE)

for (i in 1:nrow(data_score))
{
    if ((data_score[i,3]>=8)&(data_score[i,4]>=3))
    {
        text(data_score[i,1], data_score[i,2], "*", cex = 3)
    }
}

par(cex.axis = 1.5)

axis(1, 1:length(data_function_set), labels = FALSE, line = -2)

text(1:length(data_function_set)+0.4, par("usr")[3]+0.4, labels = data_function_set, srt = 45, pos = 2, xpd = TRUE, cex = 1.5)

text(length(data_function_set)/2, par("usr")[3]-11.3, labels = "Top 10% of genes             expressed in individual cell subtypes", srt = 0, pos = 1, xpd = TRUE, cex = 1.5, font = 2)

text(length(data_function_set)/2-3.1, par("usr")[3]-11.3, labels = "highly", srt = 0, pos = 1, xpd = TRUE, cex = 1.5, col = 2, font = 2)

par(cex.axis = 1.5)

axis(2, length(data_disease_set):1, data_disease_set, las = 2)

text(-11.3, 26, labels = "EGS", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

lines(c(-10.3, -10.3), c(21, 30), lty = 1, lwd = 1.5, xpd = TRUE)

text(-11, 20.5, labels = "GWAS", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

text(-11.3, 15, labels = "LC", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

lines(c(-10.3, -10.3), c(10, 19), lty = 1, lwd = 1.5, xpd = TRUE)

text(-11.3, 6.5, labels = "EGS", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

lines(c(-10.3, -10.3), c(3, 9), lty = 1, lwd = 1.5, xpd = TRUE)

text(-11.3, 2, labels = "LC", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

lines(c(-10.3, -10.3), c(1, 2), lty = 1, lwd = 1.5, xpd = TRUE)

text(-13.3, 20.5, labels = "ASD", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

lines(c(-12.3, -12.3), c(10, 30), lty = 1, lwd = 1.5, xpd = TRUE)

text(-13.3, 5.5, labels = "DD", srt = 0, pos = 1, xpd = TRUE, cex = 1.5)

lines(c(-12.3, -12.3), c(1, 9), lty = 1, lwd = 1.5, xpd = TRUE)

text((max(data_score[,1])+1)/2, par("usr")[3]+max(data_score[,2])+3.2, labels = "Telencephalon projecting excitatory neurons", srt = 0, pos = 1, xpd = TRUE, cex = 2, font = 2)

dev.off()
