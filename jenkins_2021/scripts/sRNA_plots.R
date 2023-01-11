library(tidyverse)
library(ggplot2)
library(gridExtra)

setwd("/storage/sequencing_data/paramecium_srna_sequencing/working/ben_new/scramble_map")

condition_name <- '_s3_day3'
bin_name <- 'scramble'
scaling_factor <- 12.16

# read forwards
r1_forward <- read.csv(paste('table/scramble_3062_',condition_name,'_r1_val_1.fq.gz_forward_table.txt', sep = ""), header = FALSE, sep = ",")
r2_forward <- read.csv(paste('table/scramble_3062_',condition_name,'_r2_val_2.fq.gz_forward_table.txt', sep = ""), header = FALSE, sep = ",")
# merge r1 and r2, count frequency of reads
forwards <- rbind(r1_forward, r2_forward)
forward_count <- plyr::count(forwards, 'V2')
# output dataframe to CSV
write.csv(forward_count, file = paste('forward_counts_', condition_name, '_', bin_name, ".csv", sep = ""), row.names=FALSE)
# reduce table to 20 -30 bp only
forward_count <- forward_count[-c(1,2,3,4,16,17,18,19,20,21,22,23,24,25,26,27), ]
# normalise for reads
forward_count$freq <- forward_count$freq / scaling_factor #(mean(forwards$V2))
# indicator for colour bar
forward_count <- mutate(forward_count, highlight_flag = ifelse(V2 == '23', T, F))

# read reverses
r1_reverse <- read.csv(paste('table/scramble_3062_',condition_name,'_r1_val_1.fq.gz_reverse_table.txt', sep = ""), header = FALSE, sep = ",")
r2_reverse <- read.csv(paste('table/scramble_3062_',condition_name,'_r2_val_2.fq.gz_reverse_table.txt', sep = ""), header = FALSE, sep = ",")
# merge r1 and r2, count frequency of reads
reverses <- rbind(r1_reverse, r2_reverse)
reverse_count <- plyr::count(reverses, 'V2')
# output dataframe to CSV
write.csv(reverse_count, file = paste('reverse_counts_', condition_name, '_', bin_name, ".csv", sep = ""), row.names=FALSE)
# reduce table to 20 - 30 bp only
reverse_count <- reverse_count[-c(1,2,3,4,16,17,18,19,20,21,22,23,24,25,26,27), ]
# normalise for reads
reverse_count$freq <- reverse_count$freq/ scaling_factor #(mean(reverses$V2)) 
# indicator for colour bar
reverse_count <- mutate(reverse_count, highlight_flag = ifelse(V2 == '23', T, F))

# Set scales for graphs to be properly comparable
max_forward <- max(range(forward_count$freq))
max_reverse <- max(reverse_count$freq)
y_axis_max <- 0

if ( max_forward >= max_reverse) {
  y_axis_max <- max_forward  
} else {
  y_axis_max <- max_reverse
}

# plot forward columns
fwd <- ggplot(forward_count, aes(x = V2, y= freq)) + 
  geom_col(aes(fill = highlight_flag)) + 
  xlim(19, 31) +
  ylim(0, y_axis_max) +
  scale_fill_manual(values = c('#595959', '#9999CC'), guide = FALSE) + 
  ggtitle(paste('Library ', condition_name, ' - Forward Mapped Reads, Cleaned Endosymbiont Bin')) +
  xlab("") +
  ylab("Reads per 1 Million") +
  theme_classic() +
  theme(text = element_text(size=25, colour = "gray25"),
        legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_blank(),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.title = element_text(size = 22),
        axis.line = element_line(size = 1.1, colour = "gray25"),
        axis.ticks = element_line(size = 1.1, colour = "gray25"))

# plot reverse columns
rev <- ggplot(reverse_count, aes(x = V2, y= freq)) + 
  geom_col(aes(fill = highlight_flag)) + 
  xlim(19, 31) +
  ylim(0, y_axis_max) +
  scale_fill_manual(values = c('#595959', '#9999CC'), guide = FALSE) + 
  ggtitle(paste('Library ', condition_name, ' - Reverse Mapped Reads, Cleaned Endosymbiont Bin')) +
  xlab ("Sequence Length") +
  ylab("Reads per 1 Million") +
  theme_classic() +
  theme(text = element_text(size=25, colour = "gray25"),
        legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_blank(),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.title = element_text(size = 22),
        axis.line = element_line(size = 1.1, colour = "gray25"),
        axis.ticks = element_line(size = 1.1, colour = "gray25"))

plot_out <- grid.arrange(fwd, rev, nrow = 2)

ggsave(paste(condition_name, bin_name,'plot.pdf', sep = "_"), plot_out, scale = 2, dpi = 150)