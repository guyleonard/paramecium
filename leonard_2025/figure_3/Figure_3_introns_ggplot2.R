# Load necessary libraries for data manipulation and visualization
library(Biostrings)  # For handling DNA sequences
library(dplyr)       # For data manipulation
library(ggplot2)     # For creating plots
library(ggpubr)      # For publication-ready themes
library(readr)       # For reading data
library(tidyr)       # For data tidying
library(wesanderson) # For color palettes

# Set the working directory where the files are located
#setwd("C:/Users/zool2480/Documents/Transporters/R/Introns")
setwd("~/Desktop/lab_group/paramecium/stats/introns/isoquant")

# Load DNA sequences from fasta files
sequences <- readDNAStringSet("186b_introns_for_size_freq_analysis.fas")
intron_sequences <- readDNAStringSet("intron_retentions_with_fake.fasta")

# Create data frames for intron data and intron retentions
intron_data <- data.frame(Accession = names(sequences), Length = width(sequences))
intron_retentions <- data.frame(Accession = names(intron_sequences), Length = width(intron_sequences))

# Subset data to include only introns with a length of 100 or less
new_intron_data <- intron_data %>% filter(Length <= 100)
new_intron_retentions <- intron_retentions %>% filter(Length <= 100)

# Combine intron data and intron retentions, handling many-to-many relationships
combined_data <- new_intron_data %>%
  group_by(Accession) %>%
  summarize(Length = Length[1]) %>%
  left_join(new_intron_retentions %>% group_by(Accession) %>% summarize(Length_IR = Length[1]), by = "Accession") %>%
  replace_na(list(Length_IR = 0)) # Replace NA values with 0

# Count occurrences of each intron length
length_count_table <- combined_data %>%
  count(Length) %>%
  rename_with(~ "Length_Count", .cols = "n")

# Count occurrences of each intron retention length
ir_length_count_table <- combined_data %>%
  count(Length_IR) %>%
  rename(IR_Length_Count = n, Length = Length_IR)

# Merge the length count tables and calculate the intron retention ratio
final_count_table <- length_count_table %>%
  full_join(ir_length_count_table, by = "Length") %>%
  replace_na(list(Length_Count = 0, IR_Length_Count = 0)) %>%
  filter(Length != 0) %>% # Filter out rows where Length is 0
  mutate(Ratio_IR = IR_Length_Count / Length_Count, # Calculate retention ratio
         logLengthCount = log(Length_Count)) # Calculate log of length count

# Define the Gompertz function for modeling
Gompertz <- function(x, y0, ymax, k, lag) {
  y0 + (ymax - y0) * exp(-exp(k * (lag - x) / (ymax - y0) + 1))
}

# Fit the Gompertz model to the data using non-linear least squares
Gomp1 <- nls(Ratio_IR ~ Gompertz(Length, y0, ymax, k, lag),
             data = final_count_table,
             start = list(y0 = 0, ymax = 0.4, k = 0.01, lag = 10))

# Function to predict intron retention ratios using the fitted Gompertz model
gomPredict <- function(y0, ymax, k, lag, x) {
  y0 + (ymax - y0) * exp(-exp(k * (lag - x) / (ymax - y0) + 1))
}

# Create a sequence of intron lengths for predictions
predx <- seq(0, 100, by = 1)

# Compute predictions using the fitted Gompertz model
predictionsGomp <- gomPredict(
  y0 = coef(Gomp1)[1], # Extract y0 from the fitted model
  ymax = coef(Gomp1)[2], # Extract ymax from the fitted model
  k = coef(Gomp1)[3], # Extract k from the fitted model
  lag = coef(Gomp1)[4], # Extract lag from the fitted model
  x = predx # Use the sequence of intron lengths
)

# Calculate the 95% confidence bounds for ymax
ymaxLowerBound <- coef(Gomp1)[2] - (summary(Gomp1)$coefficients[2, 2] * 1.96)
ymaxUpperBound <- coef(Gomp1)[2] + (summary(Gomp1)$coefficients[2, 2] * 1.96)

# Identify where the asymptote is reached
asympReached <- predx[which(predictionsGomp >= ymaxLowerBound)][1]

# Calculate the difference between consecutive predictions to find the inflection point
inflectionIdentify <- diff(predictionsGomp)
inflectionX <- which.max(inflectionIdentify)

# Calculate scaled log count for plotting
final_count_table <- final_count_table %>%
  mutate(scaled_count = logLengthCount / ceiling(max(logLengthCount)))

# Extract colors from the Darjeeling2 palette using wesanderson
french_dispatch_colors <- wes_palette("Darjeeling2")

# Assign colors for plotting
colour1 <- french_dispatch_colors[2] # Color for points
colour2 <- french_dispatch_colors[3] # Color for left y-axis line
colour3 <- french_dispatch_colors[5] # Color for asymptote line

# Plot intron retention ratios and Gompertz model fit using ggplot2
p <- ggplot(final_count_table, aes(x = Length)) +
  geom_point(aes(y = Ratio_IR), shape = 21, fill = colour1, color = colour1) +
  geom_line(data = data.frame(x = predx, y = predictionsGomp), aes(x = x, y = y), color = "red", size = 1) +
  geom_vline(xintercept = asympReached, linetype = "dotted", color = colour3, size = 1) +
  geom_line(aes(y = scaled_count), color = colour2, size = 1) +
  scale_y_continuous(
    name = "Intron Retention Ratio",
    position = "right",
    limits = c(0, 1),  # Set limits for the right y-axis
    breaks = seq(0, 1, by = 0.2),  # Set specific breaks for the right y-axis
    sec.axis = sec_axis(~ . * 10000, name = "Total Intron Count", 
                        labels = c(1, 10, 100, 1000, 10000))  # Manually set the breaks and labels for the left y-axis
  ) +
  scale_x_continuous(limits = c(20, 100), breaks = seq(20, 100, by = 10)) +
  labs(x = "Intron Length (bp)") +
  theme_pubr() + 
  theme(
    axis.title.y.right = element_text(color = colour1),
    axis.text.y.right = element_text(color = colour1),
    axis.line.y.right = element_line(color = colour1),
    axis.ticks.y.right = element_line(color = colour1),
    axis.title.y.left  = element_text(color = colour2),
    axis.text.y.left  = element_text(color = colour2),
    axis.line.y.left = element_line(color = colour2),
    axis.ticks.y.left = element_line(color = colour2)
  )

# Print the plot
print(p)

# Save the plot as a PDF
#ggsave("figure_3_intron_count_length_and_IR_ratio.pdf", p)

# Data frame for inflection estimation polygon
inflection_df <- data.frame(
  x = c(seq(2, 101, 1), 2),
  y = c(inflectionIdentify * (1 / max(inflectionIdentify)), 0)
)

# Data frame for ymax (asymptote) polygon
asymptote_df <- data.frame(
  x = c(0, 100, 100, 0),
  ymin = c(ymaxLowerBound, ymaxLowerBound, ymaxUpperBound, ymaxUpperBound),
  ymax = c(ymaxUpperBound, ymaxUpperBound, ymaxLowerBound, ymaxLowerBound)
)

# Data frame for predictions and confidence intervals
predictions_df <- data.frame(
  x = predx,
  y = predictionsGomp,
  ymin = rep(ymaxLowerBound, length(predx)), # Lower bound for confidence interval
  ymax = rep(ymaxUpperBound, length(predx))  # Upper bound for confidence interval
)

# Plot with additional layers for inflection and asymptote polygons using ggplot2
p2 <- ggplot(final_count_table, aes(x = Length)) +
  geom_point(aes(y = Ratio_IR), shape = 21, fill = colour1, color = "grey20", size = 3) +
  geom_line(data = predictions_df, aes(x = x, y = y), color = "red", size = 1) +
  geom_ribbon(data = predictions_df, aes(x = x, ymin = ymin, ymax = ymax),
              fill = rgb(0.8, 0.8, 1, 0.5)) +
  geom_vline(xintercept = asympReached, linetype = "dotted", color = colour3, size = 1) +
  geom_line(aes(y = scaled_count), color = colour2, size = 1) +
  geom_polygon(data = inflection_df, aes(x = x, y = y), fill = rgb(0.7, 0.7, 0.7, 0.5), color = rgb(0.7, 0.7, 0.7, 0.5)) +
  geom_vline(xintercept = inflectionX + 1, linetype = "dashed", color = "grey40", size = 1) +
  scale_y_continuous(
    name = "Intron Retention Ratio",
    position = "right",
    limits = c(0, 1),  # Set limits for the right y-axis
    breaks = seq(0, 1, by = 0.2),  # Set specific breaks for the right y-axis
    sec.axis = sec_axis(~ . * 10000, name = "Total Intron Count", 
                        labels = c(1, 10, 100, 1000, 10000))  # Manually set the breaks and labels for the left y-axis
  ) +
  scale_x_continuous(limits = c(20, 100), breaks = seq(20, 100, by = 10)) +
  labs(x = "Intron Length (bp)") +
  theme_pubr() + 
  theme(
    axis.title.y.right = element_text(color = colour1),
    axis.text.y.right = element_text(color = colour1),
    axis.line.y.right = element_line(color = colour1),
    axis.ticks.y.right = element_line(color = colour1),
    axis.title.y.left  = element_text(color = colour2),
    axis.text.y.left  = element_text(color = colour2),
    axis.line.y.left = element_line(color = colour2),
    axis.ticks.y.left = element_line(color = colour2)
  )

# Print the second plot
print(p2)

# Save the second plot as a PDF
ggsave("figure_SX_intron_count_length_and_IR_ratio_extra.pdf", p2)
