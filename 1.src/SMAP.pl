#!/usr/bin/perl

#####################################################################
# SMAP: A Perl script to merge several encapsulated postscript files
# into one single EPS file.

# Copyright (C) 2020-2021 Xusheng Wang
# Currently: xusheng.wang@und.edu

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; 

# SMAP, v 1.0.0 04/10/2021 Xusheng Wang
#####################################################################


=head1 NAME

	SMAP - A pipeline for sample matching in proteogenomics

=head1 SYNOPSIS

	perl SMAP.pl [option-1 value-1] ... [option-n value-n]

	where each option is followed by a value as follows:

	  --variant_peptide,-vf		A file containing quantitative values of variant peptides (required)
	  --genotype,-g		A genotype file used sample verification (required)
	  --output,-o			Specify an output filename	  
	  --plex			Multiplex number of the isobaric labeling approach
	  -fc				Signal to Noise ratio (optional; default is 3)
	  -nl				The upper threshold of a noise level
	  --version, -h			Print version
	  --help,-h			Print help
	  --licence,-l			Print licence
	
=head1 VERSION

	1.0.0

=cut


use warnings;
use strict;
use Getopt::Long 'HelpMessage';
use Time::Piece;
use Storable;

my $Official_Version = '1.0.0';


my $opts=GetOptions(
# Special options
  'version|v'     =>   sub { print_version(0) },  
  'help|h'     =>   sub { HelpMessage(0) },
  'licence|l'     =>   sub { print_license(0) }, 
# parameters 
  'variant_peptide|vf=s' => \my $variant_peptide_quant_file,
  'genotype|g=s' => \ my $genotype_file,
  'output|o=s' => \ my $output_file,
  'plex:i'   => \(our $plex),
  'fc:i'   => \(our $fold_change),
  'nl:i' => \ our $noise_upper,  
) or print_usage();

# die unless we got the mandatory argument
#HelpMessage(1) unless $mutant_event_file and $genotype_file;

## SMAP log file
open(LOG,">SMAP.log");

## Get current time 
my $datestring = gmtime(); 

#### Initialize parameters ##########################################
## The number multiple used in the isobaric quantification method
$plex = 11 unless ($plex);

## the minimum value of signal to noise ratio
$fold_change = 3 unless ($fold_change);

## The upper limit of the noise level
$noise_upper = 65536 unless($noise_upper);
#######################################################################

## Main program #######################################################
print "\n\tSMAP start: $datestring\n"; 
print LOG "\n\tSMAP start: $datestring\n"; 

print "\n\tReading variant peptide data\n"; 
print LOG "\n\tReading variant peptide data\n";

## Reading variant peptide data  
my ($mutant_peptide_quan_ref,$mutant_peptide_position_ref,$header_data) = Extract_variant_peptides($variant_peptide_quant_file);

print "\n\tSelecting sample-specific SNPs\n";
print LOG "\n\tSelecting sample-specific SNPs\n";
my $sample_genotype_file = "sample_specific_genotype.vcf";
select_sample_specific_genotype($genotype_file,$mutant_peptide_position_ref,$sample_genotype_file);

## Convert vcf into SNP matrix
my $SNP_matrix_file = "sample_specific_genotype.snp";
vcf2variant($sample_genotype_file,$SNP_matrix_file);

print "\n\n\tInferring sample-specific genotypes\n";
print LOG "\n\n\tInferring sample-specific genotypes\n";
my $Genotype_inferred_hash = Inferring_genotype($SNP_matrix_file,$mutant_peptide_quan_ref);

print "\n\tVerifying and correcting sample identity\n";
print LOG "\n\tVerifying and correcting sample identity\n";
genotype_matching($SNP_matrix_file,$Genotype_inferred_hash,$output_file,$header_data);

$datestring = gmtime(); 
print "\n\tSMAP end: $datestring\n"; 
print LOG "\n\tSMAP end: $datestring\n"; 

################## Subroutines ###########################################
#### read genotype data 

sub Extract_variant_peptides
{
	my $mutation_file = shift @_;
	open(MUTATION,$mutation_file) || die "cannot open the mutation file!";
	my %mutant_peptide_quan;
	my %mutant_peptide_position;
	my $header = <MUTATION>;
	chop $header;
	my @header_data=split(/\t/,$header);
	
	while(<MUTATION>) 
	{
		chomp $_;
		my @data=split(/\t/,$_);
		my $peptide = $data[0];
		my $psm = $data[2];
		
		my @quan_data = @data[4..$#data];
		my $filter = Signal_filtering(\@quan_data);
		next if($filter == 1);
		my $standardized_value = normalization(\@quan_data);
		my @merged = (@data,@$standardized_value);
		$mutant_peptide_quan{$peptide}{$psm}=join("\t",@merged);
		$mutant_peptide_position{$data[3]}=$_;
	}
	return (\%mutant_peptide_quan,\%mutant_peptide_position,\@header_data);	
}

sub Signal_filtering
{
	my $quan_data = shift;
	my $filter=0;
	my @non_zero_quan_data = @$quan_data;

	my $max_value = max(@non_zero_quan_data);
	my $min_value = min(@non_zero_quan_data);
	
	if($min_value > $noise_upper)
	{
		$filter=1;
		return $filter;
	}
#	if($max_value < $log2_min_intensity)
#	{
#		print "Lower limit filtering by max intensity: $data[0]","\r";
#		$filter=1;
#		return $filter;
#	}
	my $fold = $max_value / ($min_value+1);

	if($fold<$fold_change)
	{
		$filter=1;

		return $filter;		
	}
	return $filter;
}


## Generate sample specific genotype file

sub select_sample_specific_genotype
{
	my ($genotype_file,$mutant_peptide_position_ref,$sample_genotype_file) = @_;

	open(SAMPLE_SNP, ">$sample_genotype_file");
	open(GENO,"$genotype_file") || die "cannot open the genotype file\n";
	print "\n	Creating Batch-specific SNP matrix. It takes a while\n";

	my %SNP_mat_hash;
	my $header="";
	while(<GENO>)
	{
		chomp $_;
		if($_=~/\#CHROM/)
		{
			$header=$_;
			print SAMPLE_SNP $header,"\n";
		}
		next if($_=~/^\#\#/);
		my @data = split(/\t/,$_);
		print "\tReading SNP: $data[2]","\r";
		my $ref = $data[3];
		my $alt = $data[4];
		my $chr_pos=$data[2];
		if($mutant_peptide_position_ref->{$chr_pos})
		{
			print SAMPLE_SNP $_,"\n";
		}			
	}
	close(GENO);
	close(SAMPLE_SNP);
}

#### convert sample-specific genotype into SNP matrix 

sub vcf2variant
{
	my ($vcf_file,$snp_matrix_file) = @_;;

	open(SNPMATRIX,">$snp_matrix_file") || die "can not open the vcf file";

	
	open(VCF,"$vcf_file") || die "cannot open the input file\n";
###################################
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2014-2194
####################################
	my $header = <VCF>;
	print SNPMATRIX $header;
	
	while(<VCF>)
	{
		chomp $_;
		my @data = split(/\t/,$_);
		print SNPMATRIX join("\t",@data[0..8]);
		
		my $ref = $data[3];
		my $alt = $data[4];
		for (my $i = 9; $i<$#data;$i++)
		{
			if($data[$i] eq "0/0")
			{
				print SNPMATRIX "\t",$ref;
			}
			elsif($data[$i] eq "0/1")
			{
				print SNPMATRIX "\tH";
			}
			elsif($data[$i] eq "1/1")
			{
				print SNPMATRIX "\t",$alt;				
			}
			elsif($data[$i] eq "./.")
			{
				print SNPMATRIX "\t.";				
			}
			else
			{
				print $data[$i],"\n";
				exit;
			}
		}
		print SNPMATRIX "\n";	
	}
	close(SNPMATRIX);
}


sub Inferring_genotype
{
	my ($sample_genotype_file,$MS_Normalized_data) = @_;
	open(GENO,$sample_genotype_file);
	my %Genotype;
	my %Genotype_allele;
	while(<GENO>)
	{
		chomp $_;
		my @data = split(/\t/,$_);
		
		for(my $i=9;$i<=$#data;$i++)
		{
			next if($data[$i] eq "NA");
			next if($data[$i] eq ".");
			#my $key = $data[0] . ":" . $data[1];
			my $key = $data[2];
			$Genotype{$key}{$data[$i]}++;
			$Genotype_allele{$key}{'ref'} = $data[3];
			$Genotype_allele{$key}{'alt'} = $data[4];
		}
	}
	
	my %Genotype_hash;
			
	foreach my $peptide (keys %$MS_Normalized_data)
	{
		foreach my $PSM (keys %{$MS_Normalized_data->{$peptide}})
		{
		
			my $quan_value = $MS_Normalized_data->{$peptide}->{$PSM};
			my @data = split(/\t/,$quan_value);
			my $chr_pos = $data[3];
			my $num_genotypes = scalar keys %{$Genotype{$chr_pos}};

			next if ($num_genotypes == 0);
			next if ($num_genotypes == 1);			
		
	
			my $start = $plex + 4;
			if($num_genotypes == 2)
			{
				if(defined($Genotype{$chr_pos}{H}))
				{
					for(my $i=$start;$i<=$#data;$i++)
					{
						if($data[$i]<0.5)
						{

							my $ref_allele = $Genotype_allele{$chr_pos}{'ref'};
							$Genotype_hash{$chr_pos}{$i}{$ref_allele}++;
						}
						else
						{
							$Genotype_hash{$chr_pos}{$i}{"H"}++;
						}
					}
				}
				else
				{
					for(my $i=$start;$i<=$#data;$i++)
					{
						if($data[$i]<0.5)
						{
							my $ref_allele = $Genotype_allele{$chr_pos}{'ref'};
							$Genotype_hash{$chr_pos}{$i}{$ref_allele}++;						
						}
						else
						{
							my $alt_allele = $Genotype_allele{$chr_pos}{'alt'};
							$Genotype_hash{$chr_pos}{$i}{$alt_allele}++;						
						}
					}
				}
			}
			elsif($num_genotypes == 3)
			{
				for(my $i=$start;$i<=$#data;$i++)
				{
					if($data[$i]<0.25)
					{
						my $ref_allele = $Genotype_allele{$chr_pos}{'ref'};
						$Genotype_hash{$chr_pos}{$i}{$ref_allele}++;					
					}
					elsif($data[$i]<0.75)
					{
						$Genotype_hash{$chr_pos}{$i}{"H"}++;					
					}
					else
					{
						my $alt_allele = $Genotype_allele{$chr_pos}{'alt'};
						$Genotype_hash{$chr_pos}{$i}{$alt_allele}++;					
					}
				}
			}
			else
			{
				
				for(my $i=$start;$i<=$#data;$i++)
				{
					if($data[$i]<0.25)
					{
						my $ref_allele = $Genotype_allele{$chr_pos}{'ref'};
						$Genotype_hash{$chr_pos}{$i}{$ref_allele}++;				
					}
					elsif($data[$i]<0.75)
					{
						$Genotype_hash{$chr_pos}{$i}{"H"}++;					
					}
					else
					{
						my $alt_allele = $Genotype_allele{$chr_pos}{'alt'};
						$Genotype_hash{$chr_pos}{$i}{$alt_allele}++;					
					}
				}		
			}			
		}
	}
	my %sample_centred_inferred_hash;
	
	open(INFERGENO,">inferred_genotype.txt");
	print INFERGENO "SNP\t", join("\t",@$header_data[4..(3+$plex)]),"\n";
	foreach my $pos (keys %Genotype_hash)
	{
		print INFERGENO $pos,"\t";			
		foreach my $sample (sort {$a<=>$b} keys %{$Genotype_hash{$pos}})
		{
			foreach my $allele (reverse sort {$Genotype_hash{$pos}{$sample}{$a} <=>$Genotype_hash{$pos}{$sample}{$b}} keys %{$Genotype_hash{$pos}{$sample}})
			{
				$sample_centred_inferred_hash{$sample}{$pos}=$allele;
				print INFERGENO $allele,"\t";
				last;
			}
		}
		print INFERGENO "\n";
	}
	close(INFERGENO);
	return \%sample_centred_inferred_hash;

}

## Matching and scoring

sub genotype_matching
{
	my ($Genotype_MS_matrix,$Genotype_inferred_hash,$Final_results,$header_data) = @_;
	open(RES,">$Final_results");
	my %Genotype_MS;
	
	open(GENOTYPEMS, $Genotype_MS_matrix);
	my $header_line = <GENOTYPEMS>;
	my @header =  split(/\t/,$header_line);
	while(<GENOTYPEMS>)
	{
		next if /^\s*$/;
		chomp $_;
		
		my @data = split(/\t/,$_);
		for(my $i=9;$i<=$#data;$i++)
		{
			my $pos = $data[2];
			$Genotype_MS{$header[$i]}{$pos}=$data[$i]; 
		}
	}
		
	my %sample_match;
	my %match_ratio;
	my $i=1;
	foreach my $MS_sample (keys %$Genotype_inferred_hash)
	{
		foreach my $genotype_sample (keys %Genotype_MS)
		{
			foreach my $Genotype_inf_pos (keys %{$Genotype_inferred_hash->{$MS_sample}})
			{
				if($Genotype_inferred_hash->{$MS_sample}->{$Genotype_inf_pos} eq $Genotype_MS{$genotype_sample}{$Genotype_inf_pos})
				{
					$sample_match{$MS_sample}{$genotype_sample}{'matched'}++;
				}
				else
				{
					$sample_match{$MS_sample}{$genotype_sample}{'unmatched'}++;
				}
			}
			if($sample_match{$MS_sample}{$genotype_sample}{'unmatched'}==0)
			{
				$sample_match{$MS_sample}{$genotype_sample}{'unmatched'} = 1;
			}
			$match_ratio{$MS_sample}{$genotype_sample}=$sample_match{$MS_sample}{$genotype_sample}{'matched'} / $sample_match{$MS_sample}{$genotype_sample}{'unmatched'};
			
		}
	}
	
	print RES "Sample ID\tInferred ID\tCsore\tDeltaCscore\n";
	foreach my $MS_sample (sort {$a<=>$b} keys %match_ratio)
	{
		my $loop = 0;
		my $top_ratio = 0;
		foreach my $genotype_sample (reverse sort {$match_ratio{$MS_sample}{$a}<=>$match_ratio{$MS_sample}{$b}} keys %{$match_ratio{$MS_sample}})
		{
			$loop++;
			if($loop == 1)
			{
				$top_ratio = $match_ratio{$MS_sample}{$genotype_sample};
				print RES $$header_data[$MS_sample-$plex],"\t",$genotype_sample,"\t",$match_ratio{$MS_sample}{$genotype_sample},"\t";
			}
			if($loop == 2)
			{
				my $delta_ratio = ($top_ratio - $match_ratio{$MS_sample}{$genotype_sample}) / $top_ratio;
				print RES $delta_ratio,"\n";
				last;
			}
		}
	}
	close(RES);
	open(CSCORE,">Score.txt");
	print CSCORE "Sample ID\tInferred ID\tCsore\tDeltaCscore\n";
	foreach my $MS_sample (sort {$a<=>$b} keys %match_ratio)
	{
		my $loop = 0;
		my $top_ratio = 0;		
		foreach my $genotype_sample (reverse sort {$match_ratio{$MS_sample}{$a}<=>$match_ratio{$MS_sample}{$b}} keys %{$match_ratio{$MS_sample}})
		{
			$loop++;
			if($loop == 1)
			{
				$top_ratio = $match_ratio{$MS_sample}{$genotype_sample};
				print CSCORE $$header_data[$MS_sample-$plex],"\t",$genotype_sample,"\t",$match_ratio{$MS_sample}{$genotype_sample},"\t";
			}
			if($loop > 1)
			{
				my $delta_ratio = ($top_ratio - $match_ratio{$MS_sample}{$genotype_sample}) / $top_ratio;
				print CSCORE $delta_ratio,"\n";
				print CSCORE $$header_data[$MS_sample-$plex],"\t",$genotype_sample,"\t",$match_ratio{$MS_sample}{$genotype_sample},"\t";
			}
		}
		print CSCORE "0\n";
	}
	close(CSCORE);	
	
}


## Argument subroutines ##################################
=head
unless($variant_peptide_quant_file or $genotype_file) {
    print "	\n\tERROR: No variant peptide file or genotypic file are provided\n";
	print LOG "	\n\tERROR: No variant peptide file or genotypic file are provided\n";
    &print_usage;
  #  exit(1);
}
=cut
sub print_usage {
    print STDOUT <<USAGE;

	Usage: perl SMAP.pl [option-1 value-1] ... [option-n value-n]

	where each option is followed by a value as follows:

	  --variant_peptide,-vf		A file containing quantitative values of variant peptides (required)
	  --genotype,-g		A genotype file used sample verification (required)
	  --output,-o		Specify an output filename
	  --plex			Multiplex number of the isobaric labeling approach
	  -fc				Signal to Noise ratio (optional; default is 3)
	  -nl				The upper threshold of a noise level
	  --version,-v			Print SMAP version
	  --help,-h			Print help
	  --licence,-l		Print the licence
	For further information and options see the epsmerge documentation.
USAGE
}

sub print_version {
    print STDOUT "\n	SMAP $Official_Version\n";
}

sub print_license {
    print STDOUT <<LICENSE;

	SMAP Copyright (C) Xusheng Wang 2020-2021.
	SMAP is distributed "AS IS" under the GNU General Public License in the
	hope that it may be useful; see the file COPYING for details.
LICENSE
}

######################################################



## Basic subroutines #################################
## Standardize all signals
sub normalization
{
	my $data_ref = shift;
	
	my @data_array= @$data_ref;
	my $max = max(@data_array);
	my $min = min(@data_array);

	my @stardardized_value=();

	foreach my $data (@data_array)
	{
		my $norm = standardize($data, $min, $max);
		push(@stardardized_value,$norm);
	}
	return \@stardardized_value;
}

## Standardize a vector of values
sub standardize {
	my ($value, $min, $max) = @_; 
	return 0 if(($max - $min)==0);
	my $normalized = ($value - $min) / ($max - $min);
	return $normalized;
}

## Get the maximal value
sub max { 
 my $max = shift;
 foreach ( @_ ) { $max = $_ if $_ > $max }
 return $max;
}

## Get the minimal value
sub min { # Numbers.
 my $min = shift;
 foreach ( @_ ) { $min = $_ if $_ < $min }
 return $min;
}

## log2 transformation for all signals 
sub log2 {
   my @n = @_;
   my @t = ();
   foreach my $n (@n){
      if ($n == 0){
         $n = '0.5';
      }
      my $t = log($n)/log(2);
      $t = sprintf("%.2f",$t);
      push(@t,$t);
   }
   return(@t);
}

## calcuate the average of signals 
sub avg {
    my $total;
    $total += $_ foreach @_;
    return $total / @_;
}
######################################################