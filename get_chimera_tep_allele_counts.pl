#!/usr/bin/perl -w
##Eric Smith

use strict;
use DBI;
prepareDbh();

my $chimera_tep_counts_table = 'genomes_chimera_tep_allele_counts';

$dbh->do(qq{drop table if exists $chimera_tep_counts_table});
$dbh->do(qq{create table $chimera_tep_counts_table
            (gac_id int auto_increment primary key,
			 sample varchar(20),
			 species varchar(7),
			 sex varchar(1),
			 country varchar(20),
			 region varchar(20),
             pos int,
			 tep_allele int,
			 chimera_allele int,
			 other_allele int,
			 index(sample))});

$dbh->do(qq{drop table if exists genomes_sample_chimera_tep_info});
$dbh->do(qq{create table genomes_sample_chimera_tep_info
            (gscti_id int auto_increment primary key,
             sample varchar(20),
             total_sites int,
             hom_sites int,
             het_sites int,
             index(sample))});
			
my %snps;
my $snps_query = $dbh->prepare(qq{select pos, tep_nt, chimera_nt from chimera_tep_alns_snps});
$snps_query->execute();
while (my ($pos, $tep_nt, $chimera_nt) = $snps_query->fetchrow_array())
  {
  $snps{$pos}{tep_nt} = $tep_nt;
  $snps{$pos}{chimera_nt} = $chimera_nt;	
  }

my %sample_lookup;
my %sample_snp_info;
my $vcf_file = '/Research/gambiae/tapl/3L_1000_genomes_data/tapl_snps.vcf';
open VCF, "$vcf_file";
my $line_count = 0;
while (my $line = <VCF>)
  {
  $line_count += 1;
  if ($line_count % 100000 == 0)
    {
    my $time = localtime;
	print STDOUT "$time\t$line_count\n";	
    }
  chomp $line;

  next if $line =~ /^##/;
  if ($line =~ /^#CHROM/)
    {
    my @values = split /\t/, $line;
	my $j = 0;
	while ($j < (scalar @values))
	  {
	  my $element = $values[$j];
	  if ($element =~ /^A\D\d+/)
	    {
	    $sample_lookup[$j] = $element;
		$sample_snp_info{$element}{total_snps} = 0;
		$sample_snp_info{$element}{hom_sites} = 0;
		$sample_snp_info{$element}{het_sites} = 0;	
	    }
	  $j++;	
	  }
	}  
  else
	{
	my @values = split /\t/, $line;
	my $position = $values[1];
	foreach my $snp_pos (keys %snps)
	  {
	  if ($position == $snp_pos)
		{
		my $ref_allele = $values[3];
		my @alt_alleles = split/,/, $values[4];
		my $k = 9;
		while ($k < scalar @values)
		  {
		  my $sample_id = $sample_lookup[$k];
	      my @gt_field = split /:/, $values[$k];
		  my $genotype = $gt_field[0];
		  unless ($genotype eq './.')
			{
			$sample_snp_info{$sample_id}{total_sites} += 1;
			my @gt_split = split /\//, $genotype;
			my $allele_1 = $gt_split[0];
			my $allele_2 = $gt_split[1];
			my $tep_allele_count = 0;
			my $chimera_allele_count = 0;
			my $other_allele_count = 0;
     		my $a1_nt;
			my $a2_nt;
			if ($allele_1 == $allele_2)
			  {
			  $sample_snp_info{$sample_id}{hom_sites} += 1;	
			  }
			else
			  {
			  $sample_snp_info{$sample_id}{het_sites} += 1;	
			  }  

			if ($allele_1 == 0)
			  {
			  $a1_nt = $ref_allele;	
			  }
			else
			  {
			  $a1_nt = $alt_alleles[$allele_1 - 1];	
			  }	
  			if ($allele_2 == 0)
  			  {
  			  $a2_nt = $ref_allele;	
  			  }
  			else
  			  {
  			  $a2_nt = $alt_alleles[$allele_2 - 1];	
  			  }	
			  
			if ($a1_nt eq $snps{$snp_pos}{tep_nt})
			  {
			  $tep_allele_count += 1;
			  }
			elsif ($a1_nt eq $snps{$snp_pos}{chimera_nt})
			  {
			  $chimera_allele_count += 1;	
			  }
			else
			  {
			  $other_allele_count += 1;	
			  }		
  			if ($a2_nt eq $snps{$snp_pos}{tep_nt})
  			  {
  			  $tep_allele_count += 1;
  			  }
  			elsif ($a2_nt eq $snps{$snp_pos}{chimera_nt})
  			  {
  			  $chimera_allele_count += 1;	
  			  }
  			else
  			  {
  			  $other_allele_count += 1;	
  			  }		
				
			my ($species, $sex, $country, $region) = $dbh->selectrow_array(qq{select if(form = 'M', 'coluzzii', 'gambiae'), sex, country, region
			                                                                  from lk_genomes_project_samples
																			  where sample_id like '${sample_id}%'});
	
			$dbh->do(qq{insert into $chimera_tep_counts_table
			            values
						(null, '$sample_id', '$species', '$sex', '$country', '$region', $position, $tep_allele_count, $chimera_allele_count, $other_allele_count)});
			}	
		  $k++;
		  }	
		}	
	  }	
	}  		
  }

foreach my $sample (keys %sample_snp_info)
  {
  $dbh->do(qq{insert into genomes_sample_chimera_tep_info
              values
		      (null, '$sample', $sample_snp_info{$sample}{total_sites}, $sample_snp_info{$sample}{hom_sites}, $sample_snp_info{$sample}{het_sites})});	
  }
sub prepareDbh
  {
  my $username = #REDACTED#;
  my $password = #REDACTED#;
  my $database = #REDACTED#;
  
  
  my $dsn = "DBI:mysql:$database:localhost";
  our $dbh = DBI->connect($dsn,$username,$password);
  }
