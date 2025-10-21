# Load necessary libraries for data manipulation and visualization
library(purrr)         # For functional programming and iterating
library(ggseqlogo)     # For generating sequence logos
library(Biostrings)    # For DNA sequence manipulation

# Set the working directory to the location of input files
setwd("~/Desktop/lab_group/paramecium/stats/introns/final")

# Create a matrix containing file paths, titles, and GFF paths for sequence processing
file_title_list <- matrix(
  c(
    # Each entry contains: FASTA file path, plot title, GFF file path
    "fasta/pbursaria_186b_portal_introns.fas", "P. bursaria 186b", "gff/pbursaria_186b_portal_introns_named.gff",
    "fasta/pbursaria_110224_portal_introns.fas", "P. bursaria 110224 portal", "gff/pbursaria_110224_portal_introns_named.gff",
    "fasta/pbursaria_Dd1_portal_introns.fas", "P. bursaria DD1 portal", "gff/pbursaria_Dd1_portal_introns_named.gff",
    "fasta/pbursaria_HK1_portal_introns.fas", "P. bursaria HK1 portal", "gff/pbursaria_HK1_portal_introns_named.gff",
    "fasta/pbursaria_KM2_portal_introns.fas", "P. bursaria KM2 portal", "gff/pbursaria_KM2_portal_introns_named.gff",
    "fasta/pbursaria_STL3_portal_introns.fas", "P. bursaria STL3 portal", "gff/pbursaria_STL3_portal_introns_named.gff",
    "fasta/pcaudatum_v1_portal_introns.fas", "P. caudatum v1 portal", "gff/pcaudatum_v1_portal_introns_named.gff",
    "fasta/pcaudatum_v2_portal_introns.fas", "P. caudatum v2 portal", "gff/pcaudatum_v2_portal_introns_named.gff",
    "fasta/ptetraurelia_portal_introns.fas", "P. tetraurelia portal", "gff/ptetraurelia_portal_introns_named.gff"
    # Additional entries can be added here if needed
  ),
  ncol = 3,  # Number of columns in the matrix
  byrow = TRUE  # Fill the matrix by rows
)

# Function to create sequence names from GFF3 columns to match those in the FASTA
make_seq_name <- function(gff_row) {
  part1 <- gff_row[1]  # First part of the name from the first column
  part2 <- paste(gff_row[4] - 1, gff_row[5], sep="-")  # Second part with start and end positions
  paste(part1, part2, sep=":")  # Combine parts into a single name
}

# Function to generate sequence logos for a given FASTA file and title
generate_sequence_logos <- function(fas_file_name, title) {
  # Read intron sequences from the FASTA file
  introns <- read.fasta(fas_file_name)
  
  # Remove duplicate sequence names
  introns <- introns %>% distinct(seq.name, .keep_all = TRUE)
  
  # Calculate the length of each sequence
  introns$seq.length <- stringr::str_count(introns$seq.text)
  
  # Generate the GFF file name from the FASTA file name
  gff_file_name <- str_replace(fas_file_name, "fasta/", "gff/")
  gff_file_name <- str_replace(gff_file_name, "\\.fas$", "_named.gff")
  
  # Read the GFF3 file containing sequence annotations
  gff3 <- read.table(gff_file_name, sep="\t", comment.char="#", header=FALSE)
  
  # Remove duplicate rows based on specific columns
  gff3 <- gff3 %>% distinct(V1, V4, V5, .keep_all = TRUE)
  
  # Loop through each row in the GFF3 file to check strand and reverse complement if necessary
  for (i in 1:nrow(gff3)) {
    seq_name <- make_seq_name(gff3[i, ])  # Generate the sequence name
    intron_row <- which(introns$seq.name == seq_name)  # Find the matching intron row
    
    # If the sequence is on the negative strand, reverse complement it
    if (length(intron_row) == 1 && gff3[i, 7] == '-') {
      original_seq <- as.character(introns$seq.text[intron_row])
      complement_seq <- as.character(Biostrings::reverseComplement(DNAString(original_seq)))
      introns$seq.text[intron_row] <- complement_seq
    }
  }
  
  # Create an empty list to store filtered sequences based on length
  intron_list <- list()
  
  # Define the sequence lengths of interest (20 to 30 bp)
  seq_lengths <- 20:30
  
  # Loop through each sequence length and filter the data
  for (length in seq_lengths) {
    filtered_introns <- introns %>% filter(seq.length == length)  # Filter introns of the current length
    intron_list[[as.character(length)]] <- filtered_introns$seq.text  # Store in the list
  }
  
  # Generate the sequence logos for each filtered list
  logo_plot <- ggseqlogo(intron_list,
                         ncol = 3,  # Number of columns in the output plot
                         seq_type = "DNA",  # Type of sequence data
                         method = "prob")  # Use probability for logo
  
  # Generate the output filename by replacing non-alphanumeric characters with underscores
  output_filename <- paste0(gsub("[^[:alnum:]]", "_", title), "_logo_intron_plots.pdf")
  
  # Save the generated sequence logo plot as a PDF file
  ggsave(output_filename,
         width = 15,  # Width of the output PDF
         height = 20,  # Height of the output PDF
         dpi = 300)  # Resolution of the PDF
}

# Loop through the file_title_list and generate sequence logos for each file
walk2(file_title_list[, 1], file_title_list[, 2], generate_sequence_logos)
