args = commandArgs(TRUE)
infile <- args[1]
name_str <- args[2]
outfile <- args[3]

mytab = read.table(infile, sep = '\t', stringsAsFactors = FALSE)
valid_freqs <- mytab[, 1]

f_len <- length(valid_freqs)
umi_len <- 6
unused_len <- 4^umi_len - f_len

unused_vec <- rep(0, unused_len)

total_freq <- c(valid_freqs, unused_vec)

break_num <- max(total_freq) + 20

pdf(outfile)
hist(total_freq, break_num, xlab = 'Read count for a UMI', main = name_str)
axis(1, at=seq(0, break_num, 1))
dev.off()



