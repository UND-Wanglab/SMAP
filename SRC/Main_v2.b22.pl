#!/usr/local/bin/perl

use Storable;

## Build mutant table from VCF file
#my $vcf_input = $ARGV[0];
#vcf2mutant($vcf_input);

##change JUMPn parameter
# vi jump_g_v2.3.stage1.params
# /home/xwang4/2018/Psycho_genotype_2/Database_search/merge2.WGS-PsychChip-affy_mut.txt

##Run JUMPn
#system(qq(perl /home/yli4/development/JUMPg/v2.3.6/JUMPg-master/programs/JUMPg_v2.3.pl jump_g_v2.3.stage1.params psy_dwp_h_b16_*.raw))

##Run JUMPq
#system(qq(perl /home/yli4/development/JUMPg/v2.3.6/JUMPg-master/programs/JUMPg_v2.3.pl jump_g_v2.3.stage1.params psy_dwp_h_b16_*.raw))

## The upper limit of the minimal signals in all samples
my $log2_upper_limit_min_intensity = 65536;

## the lower limit of the maximal signals in all samples
my $log2_min_intensity = 8000;

## the lower FC value of the signals in all samples
my $FC_value = 3;


##Combine mutation with quantification

my $mutant_event_file = "/home/groupdirs/wanglab/ling/project/genotype_db_search/gnm_SNP_b22/publications/mutation_events.txt";
my $reference_peptide_bed = "/home/groupdirs/wanglab/ling/project/genotype_db_search/gnm_SNP_b22/publications/reference_peptides.bed";
my $mutant_peptide_bed = "/home/groupdirs/wanglab/ling/project/genotype_db_search/gnm_SNP_b22/publications/mutation_peptides.bed";
my $raw_quan_nonzero = "/home/groupdirs/wanglab/ling/project/genotype_db_search/quan_HH_SNP_b22/raw_quan_HH_SNP_b22_psm_nonzero.txt";
my $SNPmatrix = "/home/groupdirs/wanglab/ling/project/genotype/Input_cp/SNPmatrix.txt";
#my $MS_SNPmatrix = "/home/groupdirs/wanglab/ling/project/genotype/Input_cp/MS_SNPmatrix.txt";
my $result_file = "/home/groupdirs/wanglab/ling/project/genotype/Input_cp/Batch_ling/Batch22_Mutant_Reference_quan_results.txt";
my $Genotype_MS_matrix = "/home/groupdirs/wanglab/ling/project/genotype/Input_cp/Batch_ling/MS_b22_SNPmatrix_head.txt";
my $Final_results = "/home/groupdirs/wanglab/ling/project/genotype/Input_cp/Batch_ling/Batch22_inferred_results.txt";
print "starting extract mutant peptides\n";
my $mutant_peptide_quan_ref = Extract_mutant_peptides($mutant_event_file,$raw_quan_nonzero);
print "starting extract reference peptides\n";
$mutant_peptide_quan_ref = Extract_reference_peptides($mutant_peptide_bed,$reference_peptide_bed,$raw_quan_nonzero,$mutant_peptide_quan_ref);
print "Writing quantification results\n";
write_results($mutant_peptide_quan_ref,$result_file);
print "Extracting genotypes for this batch\n";
Extract_SNPs($SNPmatrix,$result_file,$Genotype_MS_matrix);

print "Inferring SNPs for this batch\n";
my $Genotype_inferred_hash = Inferring_genotype($Genotype_MS_matrix,$result_file);
print "Inferring samples for this batch\n";
genotype_matching($Genotype_MS_matrix,$Genotype_inferred_hash,$Final_results);

## convert vcf file containing genotypic data into mutant table for JUMPn as input
sub vcf2mutant
{
	my $vcf_file = shift;
	my $mutant_table = $vcf_file;
	$mutant_table =~s/\.vcf/\_mut\.txt/;
	open(MUT,">$mutant_table");
	print MUT "#CHROM	POS	POS	REF	ALT\n";
#####################################	
# #CHROM	POS	POS	REF	ALT
#chr1	135013	135013	A	G
#chr1	135032	135032	G	A
#chr1	135040	135040	T	C
#####################################
	
	open(VCF,"$vcf_input") || die "cannot open the input file\n";
###################################
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2014-2194
####################################

	while(<VCF>)
	{
		chomp $_;
		my @data = split(/\t/,$_);
		if($data[0] =~/^\#/)
		{
			next;
		}
		else
		{
			print MUT $data[0],"\t",$data[1],"\t",$data[1],"\t",$data[3],"\t",$data[4],"\n";
		}
	}
	return $mutant_table;
}


### Extract quantification values for mutant peptides
### First, extract it from both zero and non zero tables 
sub Extract_mutant_peptides
{
## mutant_event table from JUMPg
	my ($mutation_file,$quan_nonzero_psm_table) = @_;
	
	open(MUTATION,$mutation_file);
	my %mutant_peptide;
	my %mutant_peptide_quan;
	
	while(<MUTATION>) 
	{
		chomp $_;
		@data=split(/\t/,$_);
		$mutant_peptide{$data[4]}=$_;
	}
	my ($zero_table_ref,$nonzero_table_ref) = read_quan_data($quan_nonzero_psm_table);
	my %zero_table = %$zero_table_ref;
	my %nonzero_table = %$nonzero_table_ref;

## merging mutant peptide with quantification
	foreach my $peptide (keys %mutant_peptide)
	{
		foreach my $PSM (keys %{$zero_table{$peptide}})
		{
			$mutant_peptide_quan{$peptide}{$PSM} = $mutant_peptide{$peptide} . "\t" . $zero_table{$peptide}{$PSM} . "\tzero";
		}
		foreach my $PSM (keys %{$nonzero_table{$peptide}})
		{
			$mutant_peptide_quan{$peptide}{$PSM} = $mutant_peptide{$peptide} . "\t" . $nonzero_table{$peptide}{$PSM} . "\tnonzero";		
		}
	}
	return \%mutant_peptide_quan;
}


sub read_quan_data
{
	my ($quan_nonzero_psm_table) = shift;
	
	my $zero_table = $quan_nonzero_psm_table;
	$zero_table =~s/nonzero/zero/;
	
	open(ZERO,$zero_table);
	my %zero_table;
	<ZERO>;
	<ZERO>;
	
	while(<ZERO>)
	{
		chomp $_;
		my @data = split(/\;/,$_);
		my @peptide = split(/\./,$data[0]);
		my @quan_data = @data[27..37];
		my @non_zero_quan_data;
		foreach (@quan_data)
		{
			if($_>0)
			{			
				push(@non_zero_quan_data,$_);
			}
		}
		my $max_value = max(@non_zero_quan_data);
		my $min_value = min(@non_zero_quan_data);
		if($min_value > $log2_upper_limit_min_intensity)
		{
			print "Upper limit filtering by min intensity: $data[0]","\r";
			next;
		}
		if($max_value < $log2_min_intensity)
		{
			print "Lower limit filtering by max intensity: $data[0]","\r";
			next;
		}
		
		
		#print "Standardizing peptide in zero table: $data[0]","\r";
		my $stardardized_value = normalization(\@quan_data);		
		$zero_table{$peptide[1]}{$data[2]} = join("\t",@data[1..2]) . "\t" . join("\t",@data[27..37]) . "\t" . join("\t",@$stardardized_value);
	}
	print "\n";
	
	open(NONZERO,$quan_nonzero_psm_table);
	my %nonzero_table;
	<NONZERO>;
	<NONZERO>;
	
	while(<NONZERO>)
	{
		chomp $_;
		my @data = split(/\;/,$_);
		my @peptide = split(/\./,$data[0]);
		my @quan_data = @data[33..43];
		my @non_zero_quan_data;
		foreach (@quan_data)
		{
			if($_>0)
			{			
				push(@non_zero_quan_data,$_);
			}
		}
		my $max_value = max(@non_zero_quan_data);
		my $min_value = min(@non_zero_quan_data);
		next if($min_value > $log2_upper_limit_min_intensity);
		next if($max_value < $log2_min_intensity);		
		#print "Standardizing peptide in non_zero table: $data[0]","\r";		
		my $stardardized_value = normalization(\@quan_data);		
		$nonzero_table{$peptide[1]}{$data[2]} = join("\t",@data[1..2]) . "\t" . join("\t",@data[33..43])  . "\t" . join("\t",@$stardardized_value);
	}
	return (\%zero_table,\%nonzero_table);
}

##################################
#	1. Extract both reference peptide and mutation peptide 
#		
#	1.1 use mutant peptide to get the peptide's chromosome and start position
#		/home/xwang4/2018/Psycho_genotype/gnm_SNP_MS/publications/mutation_peptides.bed
#	
#	1.2 use the peptide's chromosome and start position to extract the corresponding reference peptides in the file of 
#		/home/xwang4/2018/Psycho_genotype/gnm_SNP_MS/publications/reference_peptides.bed
#		
#	1.3 extract quantification data from both reference and mutation peptides from the file of 
#
#		/home/xwang4/2018/Psycho_genotype/gnm_SNP_MS/intermediate/quan_Genotype/batch_1/quan_Genotype/publications/id_uni_pep_quan.txt 
##################################

sub Extract_reference_peptides
{
	my ($mutant_peptide_bed,$reference_peptide_bed,$quan_nonzero_psm_table,$mutant_peptide_quan_ref) = @_;
	
	my %reference_peptide_pos;
	my %mutant_peptide_pos;
	
	open(BED_REF,$reference_peptide_bed);
	<BED_REF>;
	while(<BED_REF>)
	{
	#chr1	100127886	100127910	K.IQEEISQK.R[1]	900	+	100127886	100127910	0	1	24,	0,
		chomp $_;
		my @data = split(/\t/,$_);
		my $peptide = $data[3];
		$peptide=~s/\[\d+\]//;
		$reference_peptide_pos{$data[0] . "_" . $data[1]} = $peptide;
	}
	close(BED_REF);	

	open(BED_MUT,$mutant_peptide_bed);
	<BED_MUT>;
	while(<BED_MUT>)
	{
	#chr1	100127886	100127910	K.IQEEISQK.R[1]	900	+	100127886	100127910	0	1	24,	0,
		chomp $_;
		my @data = split(/\t/,$_);
		my $peptide = $data[3];
		$peptide=~s/\[\d+\]//;
		$peptide=~s/(.*)\.(.*)\.(.*)/$2/;
		$mutant_peptide_pos{$peptide} = $data[0] . "_" . $data[1];
	}
	close(BED_MUT);	
	
	
	my ($zero_table_ref,$nonzero_table_ref) = read_quan_data($quan_nonzero_psm_table);
	my %zero_table = %$zero_table_ref;
	my %nonzero_table = %$nonzero_table_ref;
	
## Get corresponding 
	foreach my $peptide (keys %$mutant_peptide_quan_ref)
	{
		foreach my $mut_PSM (keys %{$mutant_peptide_quan_ref->{$peptide}})
		{	
			my @data =split(/\t/,$mutant_peptide_quan_ref->{$peptide}->{$mut_PSM});
			my $key = $mutant_peptide_pos{$peptide};
			if(defined($reference_peptide_pos{$key}))
			{
				my $reference_peptide = $reference_peptide_pos{$key};


				$reference_peptide=~s/(.*)\.(.*)\.(.*)/$2/;			
				foreach my $PSM (keys %{$zero_table{$reference_peptide}})
				{
					$mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} = $mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} . "\t" . $reference_peptide;				
					$mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} = $mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} . "\t" . $zero_table{$reference_peptide}{$PSM} . "\tzero\t";
				}
				foreach my $PSM (keys %{$nonzero_table{$reference_peptide}})
				{
					$mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} = $mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} . "\t" . $reference_peptide;								
					$mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} = $mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} . "\t" . $nonzero_table{$reference_peptide}{$PSM} . "\tnonzero\t";		
				}						
			}
			$mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} = $mutant_peptide_quan_ref->{$peptide}->{$mut_PSM} . "\n"; 
		}			
	}
	return $mutant_peptide_quan_ref;
}

sub write_results
{
	my ($mutant_peptide_quan_ref,$result_file) = @_;
	open(RESULT,">$result_file");
	
	foreach my $peptide (keys %$mutant_peptide_quan_ref)
	{
		foreach my $PSM (keys %{$mutant_peptide_quan_ref->{$peptide}})
		{	
			print RESULT $mutant_peptide_quan_ref->{$peptide}->{$PSM},"\n";
		}
	}	
}

sub normalization
{
	my $data_ref = shift;
	
	my @data_array= @$data_ref;
	my $max = max(@data_array);
	my $min = min(@data_array);

	my @stardardized_value=();

	foreach $data (@data_array)
	{
		my $norm = standardize($data, $min, $max);
		push(@stardardized_value,$norm);
	}
	return \@stardardized_value;
}

sub standardize {
	my ($value, $min, $max) = @_; 
	return 0 if(($max - $min)==0);
	$normalized = ($value - $min) / ($max - $min);
	return $normalized;
}

sub max { 
 my $max = shift;
 foreach ( @_ ) { $max = $_ if $_ > $max }
 return $max;
}

sub min { # Numbers.
 my $min = shift;
 foreach ( @_ ) { $min = $_ if $_ < $min }
 return $min;
}

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

sub avg {
    my $total;
    $total += $_ foreach @_;
    return $total / @_;
}


sub Extract_SNPs
{
	my ($SNP_matrix_file,$MS_SNP_file,$Genotype_MS_matrix) = @_;

	open(SNP, ">$Genotype_MS_matrix");
	open(SNPMAT,"$SNP_matrix_file") || die "cannot open the input file\n";
	print "Reading SNP matrix. it takes a while\n";
	my %SNP_mat_hash;
	my $header="";
	while(<SNPMAT>)
	{
		chomp $_;
		if($_=~/\#CHROM/)
		{
			$header=$_;
			last;
		}
		next if($_=~/^\#\#/);
		my @data = split(/\t/,$_);
		print "Reading SNP: $data[2]","\r";
		my $ref = $data[3];
		my $alt = $data[4];
		$SNP_mat_hash{$data[2]} = $_;

	}
	close(SNPMAT);
#	store \%SNP_mat_hash, '.SNP_mat_hash';
	
	$SNP_mat_hash_ref = retrieve('.SNP_mat_hash');
	%SNP_mat_hash = %$SNP_mat_hash_ref;

	print "\nCreating Batch-specific SNP matrix\n";
	open(MSSNP, "$MS_SNP_file");
	print SNP $header,"\n";
	while(<MSSNP>)
	{
		chomp $_;
		my @data = split(/\t/,$_);
# chr1:110709720-110709720.G_A		
		my @chr_pos_array = split(/\-/,$data[0]);
		my $chr_pos = $chr_pos_array[0];
		$chr_pos =~s/chrX/chr24/;		
		$chr_pos =~s/chr//;
		
		print SNP $SNP_mat_hash{$chr_pos},"\n";	
	}

}


sub Inferring_genotype
{
	my ($RNAseq_genotype,$MS_Normalized_data) = @_;
	open(RNASEQ,$RNAseq_genotype);
	my %RNAseq_genotype;
	while(<RNASEQ>)
	{
		chomp $_;
		my @data = split(/\t/,$_);
		
		for($i=9;$i<=$#data;$i++)
		{
			next if($data[$i] eq "NA");
			next if($data[$i] eq ".");
			my $key = $data[0] . ":" . $data[1];
			$RNAseq_genotype{$key}{$data[$i]}++;			
		}
	}
	
	my %Genotype_hash;
		
	open(MSSNP, $MS_Normalized_data);
	
	while(<MSSNP>)
	{
		next if /^\s*$/;
		chomp $_;
		
		my @data = split(/\t/,$_);
# chr1:110709720-110709720.G_A		
		my @chr_pos_array = split(/\-/,$data[0]);
		my @pos_allele = split(/\./,$chr_pos_array[1]);
		my @allele = split(/\_/,$pos_allele[1]); 
		my $chr_pos = $chr_pos_array[0];
		$chr_pos =~s/chr//;
## check how many genotype from RNAseq genotypes		
		my $num_genotypes = scalar keys %{$RNAseq_genotype{$chr_pos}};
		next if ($num_genotypes == 0);
		if(min(@data[9..19]) != 0)
		{
			next if((max(@data[9..19]) / min(@data[9..19]))< $FC_value);
		}
		print $chr_pos,"\t";		
		print join("\t",@data[1..31]),"\t";
		


		if($num_genotypes == 2)
		{
			if(defined($RNAseq_genotype{$chr_pos}{H}))
			{
				for($i=20;$i<31;$i++)
				{
					if($data[$i]<0.5)
					{
						print $allele[0],"\t";
						$Genotype_hash{$data[0]}{$i}{$allele[0]}++;
					}
					else
					{
						print "H","\t";
						$Genotype_hash{$data[0]}{$i}{"H"}++;
					}
				}
			}
			else
			{
				for($i=20;$i<31;$i++)
				{
					if($data[$i]<0.5)
					{
						print $allele[0],"\t";
						$Genotype_hash{$data[0]}{$i}{$allele[0]}++;						
					}
					else
					{
						print $allele[1],"\t";
						$Genotype_hash{$data[0]}{$i}{$allele[1]}++;						
					}
				}
			}
			print "\n";
			
		}
		elsif($num_genotypes == 3)
		{
			for($i=20;$i<31;$i++)
			{
				if($data[$i]<0.25)
				{
					print $allele[0],"\t";
					$Genotype_hash{$data[0]}{$i}{$allele[0]}++;					
				}
				elsif($data[$i]<0.75)
				{
					print "H","\t";
					$Genotype_hash{$data[0]}{$i}{"H"}++;					
				}
				else
				{
					print $allele[1],"\t";
					$Genotype_hash{$data[0]}{$i}{$allele[1]}++;					
				}
			}
			print "\n";			
		}
		else
		{
			for($i=20;$i<31;$i++)
			{
				if($data[$i]<0.25)
				{
					print $allele[0],"\t";
					$Genotype_hash{$data[0]}{$i} = $allele[0];					
				}
				elsif($data[$i]<0.75)
				{
					print "H","\t";
					$Genotype_hash{$data[0]}{$i}{"H"}++;					
				}
				else
				{
					print $allele[1],"\t";
					$Genotype_hash{$data[0]}{$i}{$allele[1]}++;					
				}
			}
			print "Not_matched\n";			
		}			
	}

	my %sample_centred_inferred_hash;
	
	foreach my $pos (keys %Genotype_hash)
	{
		print $pos,"\t";
		foreach my $sample (sort {$a<=>$b} keys %{$Genotype_hash{$pos}})
		{
			foreach $allele (reverse sort {$Genotype_hash{$pos}{$sample}{$a} <=>$Genotype_hash{$pos}{$sample}{$b}} keys %{$Genotype_hash{$pos}{$sample}})
			{
				print $allele,"\t";
				$sample_centred_inferred_hash{$sample}{$pos}=$allele;
				last;
			}
		}
		print "\n";
	}
	return \%sample_centred_inferred_hash;
}



sub genotype_matching
{
	my ($Genotype_MS_matrix,$Genotype_inferred_hash,$Final_results) = @_;
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
			my $key = $data[0] . ":" . $data[1];
			$Genotype_MS{$header[$i]}{$key}=$data[$i]; 
		}
	}
		
	my %sample_match;
	my %match_ratio;
	my $i=1;
	foreach my $MS_sample (keys %$Genotype_inferred_hash)
	{
		foreach my $genotype_sample (keys %Genotype_MS)
		{
			foreach $Genotype_inf_pos (keys %{$Genotype_inferred_hash->{$MS_sample}})
			{
				my @chr_pos_array = split(/\-/,$Genotype_inf_pos);
				my @pos_allele = split(/\./,$chr_pos_array[1]);
				my @allele = split(/\_/,$pos_allele[1]); 
				my $chr_pos = $chr_pos_array[0];
				$chr_pos =~s/chr//;				
				if($Genotype_inferred_hash->{$MS_sample}->{$Genotype_inf_pos} eq $Genotype_MS{$genotype_sample}{$chr_pos})
				{
					$sample_match{$MS_sample}{$genotype_sample}{'matched'}++;
					print $chr_pos,"\t",$MS_sample,"\t",$genotype_sample,"\t",$Genotype_inferred_hash->{$MS_sample}->{$Genotype_inf_pos},"\t",$Genotype_MS{$genotype_sample}{$chr_pos},"\t","match\n";	
				}
				else
				{
					$sample_match{$MS_sample}{$genotype_sample}{'unmatched'}++;
					print $chr_pos,"\t",$MS_sample,"\t",$genotype_sample,"\t",$Genotype_inferred_hash->{$MS_sample}->{$Genotype_inf_pos},"\t",$Genotype_MS{$genotype_sample}{$chr_pos},"\t","unmatch\n";
				}
			}
			$match_ratio{$MS_sample}{$genotype_sample}=$sample_match{$MS_sample}{$genotype_sample}{'matched'} / $sample_match{$MS_sample}{$genotype_sample}{'unmatched'};
			print RES $MS_sample,"\t",$genotype_sample,"\t",$sample_match{$MS_sample}{$genotype_sample}{'matched'},"\t",$sample_match{$MS_sample}{$genotype_sample}{'unmatched'},"\n";							
		}
	}
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
				print RES $MS_sample,"\t",$genotype_sample,"\t",$match_ratio{$MS_sample}{$genotype_sample},"\t";
			}
			if($loop == 2)
			{
				$delta_ratio = ($top_ratio - $match_ratio{$MS_sample}{$genotype_sample}) / $top_ratio;
				print RES $delta_ratio,"\n";
				last;
			}
		}
	}
	close(RES);
}
