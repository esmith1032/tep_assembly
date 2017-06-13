#!/usr/bin/perl -w
##Eric Smith

use strict;
use DBI;
prepareDbh();

my $haplotypes_table = 'genomes_samples_tep1_haplotypes';

$dbh->do(qq{drop table if exists $haplotypes_table});
$dbh->do(qq{create table $haplotypes_table
            (gac_id int auto_increment primary key,
			       sample varchar(20),
			       form varchar(1),
			       sex varchar(1),
			       country varchar(20),
			       region varchar(20),
			       ra_like int,
             rb_like int,
             r_like int,
             s_like int,
             total_sites int,
             index(sample))});

my %snp_lookup;
my $snp_query = $dbh->prepare(qq{select pos, if(ra_rb_fst > 0.5, ra_a, 0) ra_a, if(ra_rb_fst > 0.5, ra_c, 0) ra_c, if(ra_rb_fst > 0.5, ra_t, 0) ra_t, if(ra_rb_fst > 0.5, ra_g, 0) ra_g, if(ra_rb_fst > 0.5, rb_a, 0) rb_a, if(ra_rb_fst > 0.5, rb_c, 0) rb_c, if(ra_rb_fst > 0.5, rb_t, 0) rb_t, if(ra_rb_fst > 0.5, rb_g, 0) rb_g, if(ra_rb_fst < 0.5, r_a, 0) r_a, if(ra_rb_fst < 0.5, r_c, 0) r_c, if(ra_rb_fst < 0.5, r_t, 0) r_t, if(ra_rb_fst < 0.5, r_g, 0) r_g, s_a, s_c, s_t, s_g
                                 from tep1_hap_snps});
$snp_query->execute();
while (my ($pos, $ra_a, $ra_c, $ra_t, $ra_g, $rb_a, $rb_c, $rb_t, $rb_g, $r_a, $r_c, $r_t, $r_g, $s_a, $s_c, $s_t, $s_g) = $snp_query->fetchrow_array())
  {
  my %ra_hash = (A => $ra_a, C => $ra_c, T => $ra_c, G => $ra_g);  
  my %rb_hash = (A => $rb_a, C => $rb_c, T => $rb_c, G => $rb_g);  
  my %s_hash = (A => $s_a, C => $s_c, T => $s_c, G => $s_g);  
  my %r_hash = (A => $r_a, C => $r_c, T => $r_c, G => $r_g);  
  if ((sort values %ra_hash)[3] == 0)
    {
    $snp_lookup{$pos}{ra} = 'NA';  
    }
  else
    {
    $snp_lookup{$pos}{ra} = (sort {$ra_hash{$b} <=> $ra_hash{$a}} keys %ra_hash)[0]; 
    }  

  if ((sort values %rb_hash)[3] == 0)
    {
    $snp_lookup{$pos}{rb} = 'NA';  
    }
  else
    {
    $snp_lookup{$pos}{rb} = (sort {$rb_hash{$b} <=> $rb_hash{$a}} keys %rb_hash)[0]; 
    }  

  if ((sort values %s_hash)[3] == 0)
    {
    $snp_lookup{$pos}{ss} = 'NA';  
    }
  else
    {
    $snp_lookup{$pos}{ss} = (sort {$s_hash{$b} <=> $s_hash{$a}} keys %s_hash)[0]; 
    }  

  if ((sort values %r_hash)[3] == 0)
    {
    $snp_lookup{$pos}{r} = 'NA';  
    }
  else
    {
    $snp_lookup{$pos}{r} = (sort {$r_hash{$b} <=> $r_hash{$a}} keys %r_hash)[0]; 
    }  
  }

my %sample_lookup;
my %sample_allele_counts;
my $snp_count = 0;
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
  next if $line =~ /^##/;
  chomp $line;
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
		    $sample_allele_counts{$element}{r_like} = 0;
		    $sample_allele_counts{$element}{ra_like} = 0;
		    $sample_allele_counts{$element}{rb_like} = 0;
		    $sample_allele_counts{$element}{s_like} = 0;
		    $sample_allele_counts{$element}{total_sites} = 0;
        }
	    $j++;
	    }
	  }
  else
  	{
  	my @values = split /\t/, $line;
  	my $position = $values[1];
  	foreach my $snp_pos (keys %snp_lookup)
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
  			    $sample_allele_counts{$sample_id}{total_sites} += 1;
  			    my @gt_split = split /\//, $genotype;
  			    my $allele_1 = $gt_split[0];
  			    my $allele_2 = $gt_split[1];
  			    my $tep_allele_count = 0;
  			    my $chimera_allele_count = 0;
  			    my $other_allele_count = 0;
       		  my $a1_nt;
  			    my $a2_nt;
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
            
            if ($a1_nt eq $snp_lookup{$snp_pos}{ra})
  			      {
  			      $sample_allele_counts{$sample_id}{ra_like} += 1;
  			      }
  			    elsif ($a1_nt eq $snp_lookup{$snp_pos}{rb})
  			      {
    			    $sample_allele_counts{$sample_id}{rb_like} += 1;
  			      }
    			  elsif ($a1_nt eq $snp_lookup{$snp_pos}{ss})
    			    {
      			  $sample_allele_counts{$sample_id}{s_like} += 1;
    			    }
    			  elsif ($a1_nt eq $snp_lookup{$snp_pos}{r})
    			    {
      			  $sample_allele_counts{$sample_id}{r_like} += 1;
    			    }
              
            if ($a2_nt eq $snp_lookup{$snp_pos}{ra})
    			    {
    			    $sample_allele_counts{$sample_id}{ra_like} += 1;
    			    }
    			  elsif ($a2_nt eq $snp_lookup{$snp_pos}{rb})
    			    {
      			  $sample_allele_counts{$sample_id}{rb_like} += 1;
    			    }
      			elsif ($a2_nt eq $snp_lookup{$snp_pos}{ss})
      			  {
        		  $sample_allele_counts{$sample_id}{s_like} += 1;
      			  }
      			elsif ($a2_nt eq $snp_lookup{$snp_pos}{r})
      			  {
        		  $sample_allele_counts{$sample_id}{r_like} += 1;
      			  }
            }
          $k++;
          }
        }
      }
    }
  }

foreach my $sample (keys %sample_allele_counts)
  {
  my ($form, $sex, $country, $region) = $dbh->selectrow_array(qq{select form, sex, country, region
                                                                 from lk_genomes_project_samples
																                                 where sample_id like '${sample}%'});
  $dbh->do(qq{insert into $haplotypes_table
              values
			  (null, '$sample', '$form', '$sex', '$country', '$region', $sample_allele_counts{$sample}{ra_like}, $sample_allele_counts{$sample}{rb_like}, $sample_allele_counts{$sample}{r_like}, $sample_allele_counts{$sample}{s_like}, $sample_allele_counts{$sample}{total_sites})});																 
  }
  
sub prepareDbh
  {
  my $username = #REDACTED#;
  my $password = #REDACTED#;
  my $database = #REDACTED#;
  
  
  my $dsn = "DBI:mysql:$database:localhost";
  our $dbh = DBI->connect($dsn,$username,$password);
  }
