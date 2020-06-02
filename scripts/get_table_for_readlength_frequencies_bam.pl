#!/usr/bin/env perl
use warnings;
use strict;

use Bio::DB::Sam;
use Bio::Index::Fastq;
use Bio::SeqIO;
use File::Basename;
use List::MoreUtils qw(uniq);

use Data::Dumper;

my $verbose = 0;

# file names
my $sam_file   = $ARGV[0];
my $fasta_file = $ARGV[1];
my $map        = $ARGV[2];
my ( $name, $path, $suffix ) = fileparse( $fasta_file, qr/\.fasta/ );

# arrays for read matched headers
my @forward_reads_total;
my @reverse_reads_total;
my @transcripts;
my $count;

##
print "Reading: $sam_file\n";

my $bam = Bio::DB::Sam->new( -bam => "$sam_file" );
my $iterator = $bam->get_seq_stream;
while ( my $feature = $iterator->next_seq ) {
    my $read_name = $feature->display_name;
    if ( $feature->flag == 0 ) {
        push @forward_reads_total, $read_name;
    }
    elsif ( $feature->flag == 16 ) {
        push @reverse_reads_total, $read_name;
    }
}

my @forward_reads = uniq(@forward_reads_total);
my @reverse_reads = uniq(@reverse_reads_total);

#print Dumper @forward_reads;
#print Dumper @reverse_reads;

# stats
print "Done\nForward: $#forward_reads + 1 \nReverse: $#reverse_reads + 1\n";

##
print "Writing Forwards\n";

# Output lists of forward and reverse mapped reads
open my $forward_out_list, '>', "$map\_$name\_forward_list.out";
print $forward_out_list join( "\n", @forward_reads );

my $forward_fastq = `seqtk subseq $fasta_file "$map\_$name\_forward_list.out" > "$map\_$name\_forward_list.fastq"`;
my $forward_fasta = `seqtk seq -a "$map\_$name\_forward_list.fastq" > "$map\_$name\_forward_list.fasta"`;
my $forward_tidy = `rm "$map\_$name\_forward_list.fastq"`;

my @forward_matches;
my $seqio_forward =
  Bio::SeqIO->new( -file => "$map\_$name\_forward_list.fasta" );
while ( my $seq = $seqio_forward->next_seq() ) {
    my $header = $seq->id();
    my $length = $seq->length;
    push @forward_matches, "$header,$length";
}
my $forward_log = "$map\_$name\_forward_table.txt";
open my $forward_log_out, '>', $forward_log;
print $forward_log_out join( "\n", @forward_matches );
close($forward_log_out);

##
print "Writing Reverses\n";
open my $reverse_out_list, '>', "$map\_$name\_reverse_list.out";
print $reverse_out_list join( "\n", @reverse_reads );

my $reverse_fastq = `seqtk subseq $fasta_file "$map\_$name\_reverse_list.out" > "$map\_$name\_reverse_list.fastq"`;
my $reverse_fasta = `seqtk seq -a "$map\_$name\_reverse_list.fastq" > "$map\_$name\_reverse_list.fasta"`;
my $reverse_tidy = `rm "$map\_$name\_reverse_list.fastq"`;

my @reverse_matches;
my $seqio_reverse =
  Bio::SeqIO->new( -file => "$map\_$name\_reverse_list.fasta" );
while ( my $seq = $seqio_reverse->next_seq() ) {
    my $header = $seq->id();
    my $length = $seq->length;
    push @reverse_matches, "$header,$length";
}
my $reverse_log = "$map\_$name\_reverse_table.txt";
open my $reverse_log_out, '>', $reverse_log;
print $reverse_log_out join( "\n", @reverse_matches );
close($reverse_log_out);
