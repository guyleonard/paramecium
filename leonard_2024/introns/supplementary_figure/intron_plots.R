# Load the required libraries for data processing and visualization
library(phylotools)    # For handling sequence data
library(stringr)       # For string manipulation
library(dplyr)         # For data manipulation
library(ggplot2)       # For plotting
library(cowplot)       # For combining multiple plots
library(wesanderson)   # For color palettes
library(patchwork)     # For advanced plot layout management

# Function to generate and save the plot
generate_plot <- function(filename, title_name) {
  # Set the working directory to where intron data files are located
  setwd("~/Desktop/lab_group/paramecium/stats/introns")
  
  # Read the intron data from the FASTA file
  introns <- read.fasta(filename)
  
  # Calculate the length of each intron sequence
  introns$seq.length <- str_count(introns$seq.text)
  
  # Calculate the mean and standard deviation of the intron lengths
  mean_value <- round(mean(introns$seq.length), 2)
  std_dev <- round(sd(introns$seq.length), 2)
  
  # Prepare a text label for the mean and standard deviation
  stats_label_text <- substitute(
    paste("Mean: ", mean_value, " ", sigma, " ", std_dev),
    list(mean_value = mean_value, std_dev = std_dev)
  )
  # Here, sigma is a special symbol for "σ" in the expression
  
  # Group introns by length and count the number of occurrences for each length
  introns_count <- introns %>%
    group_by(seq.length) %>%
    tally(name = "count")
  
  # Filter introns with lengths between 20 and 30 for an inset plot
  introns_count_20to30 <- introns_count %>%
    filter(seq.length >= 20 & seq.length <= 30)
  
  # Create a dataset for data labels at specific lengths
  data_labels <- introns_count %>%
    filter(seq.length %in% c(20:30) |
             (seq.length > 30 & seq.length %% 10 == 0))
  
  # Find the intron length with the highest count for highlighting
  highlight_length <- introns_count$seq.length[which.max(introns_count$count)]
  
  # Create the main plot with a histogram of intron lengths
  main_plot <- ggplot() +
    # Add bars for all intron lengths
    geom_bar(
      data = introns_count,
      aes(x = seq.length, y = count),
      stat = "identity",
      fill = wes_palette("Darjeeling2", 3)[1],  # Set fill color
      alpha = 0.7  # Set transparency
    ) +
    # Highlight the bar with the highest count
    geom_bar(
      data = filter(introns_count, seq.length == highlight_length),
      aes(x = seq.length, y = count),
      stat = "identity",
      fill = wes_palette("Darjeeling2", 3)[2],  # Different color for highlight
      alpha = 0.7
    ) +
    # Add a line to connect the top of the bars
    geom_line(
      data = introns_count,
      aes(x = seq.length, y = count),
      color = wes_palette("Darjeeling2", 3)[3],
      size = 1  # Line thickness
    ) +
    # Add text labels to selected bars
    geom_text(
      data = data_labels,
      aes(x = seq.length, y = count, label = count),
      vjust = ifelse(data_labels$seq.length > 30, -0.5, 1.5),
      size = 3,  # Text size
      angle = ifelse(data_labels$seq.length > 30, 45, 0)  # Text angle
    ) +
    # Add a vertical line at the mean length
    geom_segment(
      aes(
        x = mean_value,
        y = 0,
        xend = mean_value,
        yend = min(1000, max(introns_count$count))
      ),
      color = wes_palette("Moonrise2", 4)[1],
      linetype = "solid",
      size = 1
    ) +
    # Annotate the plot with the mean and standard deviation
    annotate(
      "text",
      x = mean_value + 2,
      y = min(1000, max(introns_count$count)) / 2,
      label = stats_label_text,
      hjust = -0.05,
      vjust = 0.5,
      color = wes_palette("Moonrise2", 4)[1],
      size = 3
    ) +
    # Add labels and titles
    labs(
      x = "Intron Length (bp)",
      y = "Count",
      title = paste("Histogram of", title_name, "Intron Lengths")
    ) +
    # Customize x-axis breaks and limits
    scale_x_continuous(limits = c(19, 100),
                       breaks = seq(20, 100, by = 10)) +
    # Use a minimal theme for the plot
    theme_minimal() +
    theme(text = element_text(size = 12),
          panel.grid = element_blank()) +
    # Adjust plot clipping
    coord_cartesian(clip = "off")
  
  # Create an inset plot for lengths 20 to 30
  inset_plot <- ggplot(data = introns_count_20to30) +
    # Add bars for the inset data
    geom_bar(
      aes(x = seq.length, y = count),
      stat = "identity",
      fill = wes_palette("Darjeeling2", 3)[1],
      alpha = 0.7
    ) +
    # Highlight the bar with the highest count
    geom_bar(
      data = filter(introns_count_20to30, seq.length == highlight_length),
      aes(x = seq.length, y = count),
      stat = "identity",
      fill = wes_palette("Darjeeling2", 3)[2],
      alpha = 0.7
    ) +
    # Add text labels to the inset bars
    geom_text(aes(x = seq.length, y = count, label = count),
              vjust = -0.5,
              size = 3) +
    # Customize axis labels
    labs(x = "Intron Length", y = "Count", title = NULL) +
    # Customize x-axis breaks for the inset
    scale_x_continuous(breaks = seq(20, 30, by = 1)) +
    # Use a minimal theme for the inset plot
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10)
    )
  
  # Combine main plot and inset plot
  plot_combined <- ggdraw() +
    draw_plot(main_plot) +  # Draw the main plot
    draw_plot(
      inset_plot,  # Draw the inset plot
      x = 0.35,  # Position of the inset plot
      y = 0.25,  # Position of the inset plot
      width = 0.5,  # Width of the inset plot
      height = 0.5  # Height of the inset plot
    )
  
  # Save the combined plot as a PDF file with the title name
  output_filename <- gsub(" ", "", title_name)
  ggsave(
    paste0(output_filename, ".pdf"),  # Output file name
    plot = plot_combined,  # Plot to save
    width = 10,  # Width of the output PDF
    height = 8  # Height of the output PDF
  )
  
  # Return the combined plot
  return(plot_combined)
}

# List of filenames and title names for generating plots
file_title_list <- list(
  c("186b_introns.fas", "P. bursaria 186b Final Gene Predictions"),
  c("110224_introns.fas", "P. bursaria 110224 GeneMark ES"),
  c("110224_portal_introns.fas", "P. bursaria 110224 Genome Portal"),
  c("DD1_introns.fas", "P. bursaria DD1 GeneMark ES"),
  c("pbursaria_Dd1_portal_introns.fas", "P. bursaria DD1 Genome Portal"),
  c("HK1_introns.fas", "P. bursaria HK1 GeneMark ES"),
  c("pbursaria_HK1_portal_introns.fas", "P. bursaria HK1 Genome Portal"),
  c("KM2_introns.fas", "P. bursaria KM2 GeneMark ES"),
  c("pbursaria_KM2_portal_introns.fas", "P. bursaria KM2 Genome Portal"),
  c("STL3_introns.fas", "P. bursaria STL3 GeneMark ES"),
  c("pbursaria_STL3_portal_introns.fas", "P. bursaria STL3 Genome Portal"),
  c("caudatum_introns.fas", "P. caudatum GeneMark ES"),
  c("caudatum_v1_portal_introns.fas", "P. caudatum v1 Genome Portal"),
  c("caudatum_v2_portal_introns.fas", "P. caudatum v2 Genome Portal"),
  c("tetraurelia_introns.fas", "P. tetraurelia GeneMark ES"),
  c("ptetraurelia_portal_introns.fas", "P. tetraurelia Genome Portal")
)

# Generate individual plots and save them as PDFs with consistent y-axis limits
individual_plots <- lapply(file_title_list, function(x) {
  plot <- generate_plot(x[1], x[2])
  # Set uniform limits for the y-axis across all plots
  plot + scale_y_continuous(limits = c(0, 150))
})

# Combine all individual plots into one large plot with consistent y-axis limits
combined_plot <- cowplot::plot_grid(plotlist = individual_plots, ncol = 1)

# Save the combined plot as a single PDF file
ggsave(
  "combined_plots.pdf",  # Output file name
  plot = combined_plot,  # Combined plot to save
  width = 10,  # Width of the output PDF
  height = 50,  # Height of the output PDF
  limitsize = FALSE  # Disable size limits for the output
)
