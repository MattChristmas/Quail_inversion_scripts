#!/usr/bin/perl -w


# Script for calculating dxy outside inverted region between quails with and without the inversion.

use strict; #enforce some good programming rules
use warnings;
use diagnostics;
use Data::Dumper; 


open IN, "quail_no_inversion.recode.vcf" or die;

open OUT, ">>quail_outside_inversion_dxy.txt" or die;

print OUT "inverted_sample\tnon_inverted_sample\tDXY\n";


my %samples = ( # Set up a sample array where the key is the position of the sample in the vcf file
	9 => '83PSIII',
	10 => '78PSIII',
	11 => '67PSIII',
	12 => '65PSIII',
	13 => '8PSII',
	14 => '5PSb',
	15 => '52PSII',
	16 => '21PS',
	17 => '46PS',
	18 => '24PS',
	19 => '28PSIII',
	20 => '41PSb',
	21 => '27PSIII',
	22 => '33PSIII',
	23 => '57PSIII',
	24 => '40PS',
	);

my %inverted = (
	'8PSII' => 0,
	'5PSb' => 0,
	'52PSII' => 0,
	'21PS' => 0,
	'46PS' => 0,
	'28PSIII' => 0,
	'41PSb' => 0,
	'27PSIII' => 0,
	'33PSIII' => 0,
	'57PSIII' => 0,
	);

my %non_inverted = (
	'83PSIII' => 0,
	'78PSIII' => 0,
	'67PSIII' => 0,
	'65PSIII' => 0,
	'24PS' => 0,
	'40PS' => 0,
	);


## Keep count of differences between each haplotype comparison
my %dxy;

# Ungapped length of genome = 917,263,224
# Inversion length = 116510000
my $genome_length = 800753224; # Length of the ungapped genome minus the inversion


while (my $line = <IN>) {
	
	chomp $line;				

	if ( $line =~ m/^\##/ ) {				
		next;
	}
	
	if ( $line =~ m/^\#CHROM/ ) { 
		
		my @headerLine = split(/\t/,$line);				
	
		next;
		
	}
	
	my @elements = split(/\t/,$line); 					


	my $ref_base = $elements[3];
	my $alt_base = $elements[4];
	
	for (my $i=9; $i <= 24; $i++) { # loop through all (i) of the samples in the vcf...
		my $sampleA = $samples{$i};
		
		for (my $j=9; $j <= 24; $j++) { # and compare them to all other (j) samples
			my $sampleB = $samples{$j};
				
			if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
				next;
			}
				
			if ((exists $inverted{$sampleA}) && (exists $non_inverted{$sampleB})) { # To ensure we are only comparing samples between the two groups
				
				my $comp = "$sampleA\t$sampleB";
			
				my $data_i = $elements[ $i ];		
				my @i_data = split(/:/, $data_i);
				my $geno_i = $i_data[ 0 ]; # Get the genotype data for the ith sample
	
				my $data_j = $elements[ $j ];
				my @j_data = split(/:/, $data_j);
				my $geno_j = $j_data[ 0 ]; # get the genotype data for the jth sample
			
				if (($geno_i eq ".") || ($geno_j eq ".")) {
			
					next; # skip missing data
				}
			
				else {
			
					my @i_alleles = split(/\//,$geno_i); # Get the two alleles for the ith sample from the genotype field
					my @i_alleles_sorted = sort (@i_alleles); # This will ensure we are always comparing alleles in alphabetical order
					my $i_a1 = $i_alleles_sorted[0];
					my $i_a2 = $i_alleles_sorted[1];

	
					my @j_alleles = split(/\//,$geno_j);
					my @j_alleles_sorted = sort (@j_alleles);
					my $j_a1 = $j_alleles_sorted[0];
					my $j_a2 = $j_alleles_sorted[1];
	
					if ($i_a1 != $j_a1) { # Check whether the alleles match and if not add a count to the $hap1hap1 counter
						$dxy{$comp}++;
					}
	
					if ($i_a2 != $j_a2) {
						$dxy{$comp}++;
					}
				}	
			}
				
			else {
				next;
			}
		}	
	}
		
}
	

my $total_dxy = 0;
my $total_comparisons = 0;
	

while( my( $key, $value ) = each %dxy ){
   
   	my $dxy = $value/$genome_length;
   
   	my @samples = split(/\t/, $key);
  	my $name1 = $samples[0];
   	my $name2 = $samples[1];
   	
   	if (exists $inverted{$name1}) {
   	
   		print OUT "$name1\t$name2\t$dxy\n";
   	
   	
   	}
   	
   	else {
   	
   		print OUT "$name2\t$name1\t$dxy\n";
   	
   	}
   				
   	
   
   	$total_comparisons++;
   	$total_dxy = $total_dxy + $dxy;
	
}


my $average_dxy = $total_dxy/$total_comparisons;

print "Average dxy between inverted and non-inverted samples across genome = $average_dxy\n";
