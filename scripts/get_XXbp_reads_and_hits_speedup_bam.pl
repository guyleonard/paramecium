#!/usr/bin/env perl
use warnings;
use strict;

use Bio::SeqIO;
use Bio::DB::Sam;
use File::Basename;

use Data::Dumper;

# file names
my $fasta_file   = $ARGV[0];
my $sam_file     = $ARGV[1];
my $transcripts  = $ARGV[2];
my $input_length = $ARGV[3];
my ( $name, $path, $suffix ) = fileparse( $fasta_file, qr/\.fasta/ );

#
my @twenty3_bp_list;

## FORWARD
#fasta_file= endo_2855_n_3_r1_val_1_forward_table.txt
print "Reading: $fasta_file\t";
my $seqio_forward_r1 = Bio::SeqIO->new( -file => "$fasta_file" );
while ( my $seq = $seqio_forward_r1->next_seq() ) {
    my $header = $seq->id();
    my $length = $seq->length;
    if ( $length eq $input_length ) {
        push @twenty3_bp_list, "$header";
    }
}
print "done.\n";

$fasta_file =~ s/r1_/r2_/g;
$fasta_file =~ s/val_1/val_2/g;
print "Reading: $fasta_file\t";
#fasta_file= endo_2855_n_3_r2_val_2_forward_table.txt
my $seqio_forward_r2 = Bio::SeqIO->new( -file => "$fasta_file" );
while ( my $seq = $seqio_forward_r2->next_seq() ) {
    my $header = $seq->id();
    my $length = $seq->length;
    if ( $length eq $input_length ) {
        push @twenty3_bp_list, "$header";
    }
}
print "done.\n";

## REVERSE
$fasta_file =~ s/forward/reverse/g;
print "Reading: $fasta_file\t";

#fasta_file= endo_2855_n_3_r2_val_2_reverse_table.txt
my $seqio_reverse_r2 = Bio::SeqIO->new( -file => "$fasta_file" );
while ( my $seq = $seqio_reverse_r2->next_seq() ) {
    my $header = $seq->id();
    my $length = $seq->length;
    if ( $length eq $input_length ) {
        push @twenty3_bp_list, "$header";
    }
}
print "done.\n";

$fasta_file =~ s/r2_/r1_/g;
$fasta_file =~ s/val_2/val_1/g;
print "Reading: $fasta_file\t";
#fasta_file= endo_2855_n_3_r1_val_1_reverse_table.txt
my $seqio_reverse_r1 = Bio::SeqIO->new( -file => "$fasta_file" );
while ( my $seq = $seqio_reverse_r1->next_seq() ) {
    my $header = $seq->id();
    my $length = $seq->length;
    if ( $length eq $input_length ) {
        push @twenty3_bp_list, "$header";
    }
}
print "done.\n";

print "Reading: $sam_file\n";
my %transcript_counts;
my $bam = Bio::DB::Sam->new( -bam => "$sam_file" );
my $iterator = $bam->get_seq_stream;
while ( my $feature = $iterator->next_seq ) {
    my $read_name = $feature->display_name;
    my $transcript_ID = $feature->seq_id;
    #print "$read_name & $transcript_ID\n"; 
    if ( $feature->flag == 0 ) {
        if ( grep /$read_name/, @twenty3_bp_list ) {
		$transcript_counts{$transcript_ID}++;
	}
    }
    elsif ( $feature->flag == 16 ) {
	if ( grep /$read_name/, @twenty3_bp_list ) {
		$transcript_counts{$transcript_ID}++;
	}
    }
}


$sam_file =~ s/r1_/r2_/g;
$sam_file =~ s/val_1/val_2/g;
print "Reading: $sam_file\n";
$bam = Bio::DB::Sam->new( -bam => "$sam_file" );
$iterator = $bam->get_seq_stream;
while ( my $feature = $iterator->next_seq ) {
    my $read_name = $feature->display_name;
    my $transcript_ID = $feature->seq_id;
    if ( $feature->flag == 0 ) {
        if ( grep /$read_name/, @twenty3_bp_list ) {
                $transcript_counts{$transcript_ID}++;
        }
    }
    elsif ( $feature->flag == 16 ) {
        if ( grep /$read_name/, @twenty3_bp_list ) {
                $transcript_counts{$transcript_ID}++;
        }
    }
}

print "Output: Transcript ID List - $name\_$input_length\.list\n";
open my $list_out, '>>', "$name\_$input_length\.list";
foreach (sort keys %transcript_counts) {
    print $list_out "$_\n";
}
close $list_out;

print "Output: Transcript ID Count List - $name\_$input_length\.counts\n";
open my $count_out, '>>', "$name\_$input_length\.counts";
foreach (sort keys %transcript_counts) {
    print $count_out "$_\t$transcript_counts{$_}\n";
}
close $count_out;

print "Output: Transcript ID FASTAs $name\_$input_length\_transcripts.fasta\n";
my $transcript_fasta =
`faSomeRecords $transcripts "$name\_$input_length\.list" "$name\_$input_length\_transcripts.fasta"`;
