#!/usr/bin/perl -w
##Eric Smith

use strict;
use DBI;
prepareDbh();

my $allele_counts_table = 'genomes_allele_counts';

$dbh->do(qq{drop table if exists $allele_counts_table});
$dbh->do(qq{create table $allele_counts_table
            (gac_id int auto_increment primary key,
			 sample varchar(20),
			 form varchar(1),
			 sex varchar(1),
			 country varchar(20),
			 region varchar(20),
			 r_like int,
			 s_like int,
			 n_like int,
			 rs_like int,
			 rn_like int,
			 sn_like int,
             r_het int,
			 s_het int,
			 n_het int,
			 rs_het int,
			 rn_het int,
			 sn_het int,
			 new int,
			 index(sample))});

my %sample_lookup;
my %sample_allele_counts;
my $snp_count = 0;
my $vcf_file = '3L_1000_genomes_data/ag1000g.phase1.AR2.3L.vcf';
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
		$sample_allele_counts{$element}{s_like} = 0;
		$sample_allele_counts{$element}{n_like} = 0;
		$sample_allele_counts{$element}{rs_like} = 0;
		$sample_allele_counts{$element}{rn_like} = 0;
		$sample_allele_counts{$element}{sn_like} = 0;
		$sample_allele_counts{$element}{unique} = 0;
		$sample_allele_counts{$element}{r_het} = 0;
		$sample_allele_counts{$element}{s_het} = 0;
		$sample_allele_counts{$element}{n_het} = 0;
		$sample_allele_counts{$element}{rs_het} = 0;
		$sample_allele_counts{$element}{rn_het} = 0;
		$sample_allele_counts{$element}{sn_het} = 0;
		}
	  $j++;
	  }
	}
  else
    {
	my @values = split /\t/, $line;
	my $position = $values[1];
	
	my $alleles_query = $dbh->prepare(qq{select * from new_allele_snps where position = $position});
	$alleles_query->execute();
	my ($nas_id, $no_pos, $r_allele, $s_allele, $new_allele) = $alleles_query->fetchrow_array();
	
	if ($nas_id)
	  {
	  $snp_count += 1;
	  my $ref_allele = $values[3];
	  my @alt_alleles = split /,/, $values[4];
	  my $k = 9;
	  while ($k < scalar @values)
	    {
	    my $sample_id = $sample_lookup[$k];
	    my $sample_allele;
		my @gt_field = split /:/, $values[$k];
		my $genotype = $gt_field[0];
	    if ($genotype eq './.' or $genotype eq '0/0')
		  {
		  $sample_allele = $ref_allele;
		  if ($sample_allele eq $r_allele)
		    {
			if ($sample_allele eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_like} += 1;			  
			  }
			elsif($sample_allele eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_like} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_like} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
		    if ($sample_allele eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_like} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_like} += 1;	
			  }  	
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		    $sample_allele_counts{$sample_id}{n_like} += 1;	
		    }
		  else
		    {
		    $sample_allele_counts{$sample_id}{unique} += 1;	
		    }			
		  }
		elsif ($genotype eq '1/1')
		  {
		  $sample_allele = $alt_alleles[0];	
		  if ($sample_allele eq $r_allele)
		    {
			if ($sample_allele eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_like} += 1;			  
			  }
			elsif($sample_allele eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_like} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_like} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
		    if ($sample_allele eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_like} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_like} += 1;	
			  }  	
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		    $sample_allele_counts{$sample_id}{n_like} += 1;	
		    }
		  else
		    {
		    $sample_allele_counts{$sample_id}{unique} += 1;	
		    }					  
		  }  
  		elsif ($genotype eq '2/2')
  		  {
  		  $sample_allele = $alt_alleles[1];	
  		  if ($sample_allele eq $r_allele)
  		    {
  			if ($sample_allele eq $s_allele)
  			  {
  			  $sample_allele_counts{$sample_id}{rs_like} += 1;			  
  			  }
  			elsif($sample_allele eq $new_allele)
  			  {
  			  $sample_allele_counts{$sample_id}{rn_like} += 1;	
  			  }
  			else
  			  {
  			  $sample_allele_counts{$sample_id}{r_like} += 1;	
  			  }    
  			}
  		  elsif ($sample_allele eq $s_allele)
  		    {
  		    if ($sample_allele eq $new_allele)
  			  {
  			  $sample_allele_counts{$sample_id}{sn_like} += 1;	
  			  }
  			else
  			  {
  			  $sample_allele_counts{$sample_id}{s_like} += 1;	
  			  }  	
  		    }
  		  elsif ($sample_allele eq $new_allele)
  		    {
  		    $sample_allele_counts{$sample_id}{n_like} += 1;	
  		    }
  		  else
  		    {
  		    $sample_allele_counts{$sample_id}{unique} += 1;	
  		    }					  
  		  }
  		elsif ($genotype eq '3/3')
  		  {
  		  $sample_allele = $alt_alleles[2];	
  		  if ($sample_allele eq $r_allele)
  		    {
  			if ($sample_allele eq $s_allele)
  			  {
  			  $sample_allele_counts{$sample_id}{rs_like} += 1;			  
  			  }
  			elsif($sample_allele eq $new_allele)
  			  {
  			  $sample_allele_counts{$sample_id}{rn_like} += 1;	
  			  }
  			else
  			  {
  			  $sample_allele_counts{$sample_id}{r_like} += 1;	
  			  }    
  			}
  		  elsif ($sample_allele eq $s_allele)
  		    {
  		    if ($sample_allele eq $new_allele)
  			  {
  			  $sample_allele_counts{$sample_id}{sn_like} += 1;	
  			  }
  			else
  			  {
  			  $sample_allele_counts{$sample_id}{s_like} += 1;	
  			  }  	
  		    }
  		  elsif ($sample_allele eq $new_allele)
  		    {
  		    $sample_allele_counts{$sample_id}{n_like} += 1;	
  		    }
  		  else
  		    {
  		    $sample_allele_counts{$sample_id}{unique} += 1;	
  		    }					  
  		  }
		elsif ($genotype eq '0/1')
		  {
		  $sample_allele = $ref_allele;	
		  my $allele_2 = $alt_alleles[0];
		  if ($sample_allele eq $r_allele)
		    {
			if ($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }    
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		  	if ($allele_2 eq $r_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{rn_het} += 1;			  
		  	  }
		  	elsif($allele_2 eq $s_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{sn_het} += 1;	
		  	  }
		  	else
		  	  {
		  	  $sample_allele_counts{$sample_id}{n_het} += 1;	
		  	  }    
		    }
		  else
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;			  
			  }
			elsif($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }
			elsif ($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{n_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{unique} += 1;	
			  }      
		    }
		  }	 
		elsif ($genotype eq '0/2')
		  {
		  $sample_allele = $ref_allele;	
		  my $allele_2 = $alt_alleles[1];
		  if ($sample_allele eq $r_allele)
		    {
			if ($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }    
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		  	if ($allele_2 eq $r_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{rn_het} += 1;			  
		  	  }
		  	elsif($allele_2 eq $s_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{sn_het} += 1;	
		  	  }
		  	else
		  	  {
		  	  $sample_allele_counts{$sample_id}{n_het} += 1;	
		  	  }    
		    }
		  else
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;			  
			  }
			elsif($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }
			elsif ($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{n_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{unique} += 1;	
			  }      
		    } 
		  }	
		elsif ($genotype eq '0/3')
		  {
		  $sample_allele = $ref_allele;	
		  my $allele_2 = $alt_alleles[2];
		  if ($sample_allele eq $r_allele)
		    {
			if ($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }    
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		  	if ($allele_2 eq $r_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{rn_het} += 1;			  
		  	  }
		  	elsif($allele_2 eq $s_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{sn_het} += 1;	
		  	  }
		  	else
		  	  {
		  	  $sample_allele_counts{$sample_id}{n_het} += 1;	
		  	  }    
		    }
		  else
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;			  
			  }
			elsif($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }
			elsif ($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{n_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{unique} += 1;	
			  }      
		    }
		  }	 
		elsif ($genotype eq '1/2')
		  {
		  $sample_allele = $alt_alleles[0];	
		  my $allele_2 = $alt_alleles[1];
		  if ($sample_allele eq $r_allele)
		    {
			if ($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }    
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		  	if ($allele_2 eq $r_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{rn_het} += 1;			  
		  	  }
		  	elsif($allele_2 eq $s_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{sn_het} += 1;	
		  	  }
		  	else
		  	  {
		  	  $sample_allele_counts{$sample_id}{n_het} += 1;	
		  	  }    
		    }
		  else
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;			  
			  }
			elsif($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }
			elsif ($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{n_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{unique} += 1;	
			  }      
		    }
		  }	 
		elsif ($genotype eq '1/3')
		  {
		  $sample_allele = $alt_alleles[0];	
		  my $allele_2 = $alt_alleles[2];
		  if ($sample_allele eq $r_allele)
		    {
			if ($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }    
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		  	if ($allele_2 eq $r_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{rn_het} += 1;			  
		  	  }
		  	elsif($allele_2 eq $s_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{sn_het} += 1;	
		  	  }
		  	else
		  	  {
		  	  $sample_allele_counts{$sample_id}{n_het} += 1;	
		  	  }    
		    }
		  else
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;			  
			  }
			elsif($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }
			elsif ($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{n_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{unique} += 1;	
			  }      
		    }
		  }	 
		elsif ($genotype eq '2/3')
		  {
		  $sample_allele = $alt_alleles[1];	
		  my $allele_2 = $alt_alleles[2];
		  if ($sample_allele eq $r_allele)
		    {
			if ($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{rn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;	
			  }    
			}
		  elsif ($sample_allele eq $s_allele)
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{rs_het} += 1;			  
			  }
			elsif($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{sn_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }    
		    }
		  elsif ($sample_allele eq $new_allele)
		    {
		  	if ($allele_2 eq $r_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{rn_het} += 1;			  
		  	  }
		  	elsif($allele_2 eq $s_allele)
		  	  {
		  	  $sample_allele_counts{$sample_id}{sn_het} += 1;	
		  	  }
		  	else
		  	  {
		  	  $sample_allele_counts{$sample_id}{n_het} += 1;	
		  	  }    
		    }
		  else
		    {
			if ($allele_2 eq $r_allele)
			  {
			  $sample_allele_counts{$sample_id}{r_het} += 1;			  
			  }
			elsif($allele_2 eq $s_allele)
			  {
			  $sample_allele_counts{$sample_id}{s_het} += 1;	
			  }
			elsif ($allele_2 eq $new_allele)
			  {
			  $sample_allele_counts{$sample_id}{n_het} += 1;	
			  }
			else
			  {
			  $sample_allele_counts{$sample_id}{unique} += 1;	
			  }      
		    }
		  }	 		
		$k++;
		}
	  }	
	} 	
  }
print STDOUT "$snp_count\n";
foreach my $sample (keys %sample_allele_counts)
  {
  my ($form, $sex, $country, $region) = $dbh->selectrow_array(qq{select form, sex, country, region
                                                                 from lk_genomes_project_samples
																 where sample_id like '${sample}%'});
  $dbh->do(qq{insert into $allele_counts_table
              values
			  (null, '$sample', '$form', '$sex', '$country', '$region', $sample_allele_counts{$sample}{r_like}, $sample_allele_counts{$sample}{s_like}, $sample_allele_counts{$sample}{n_like}, $sample_allele_counts{$sample}{rs_like}, $sample_allele_counts{$sample}{rn_like}, $sample_allele_counts{$sample}{sn_like}, $sample_allele_counts{$sample}{r_het}, $sample_allele_counts{$sample}{s_het}, $sample_allele_counts{$sample}{n_het}, $sample_allele_counts{$sample}{rs_het}, $sample_allele_counts{$sample}{rn_het}, $sample_allele_counts{$sample}{sn_het}, $sample_allele_counts{$sample}{unique})});																 
  }
  
sub prepareDbh
  {
  my $username = #REDACTED#;
  my $password = #REDACTED#;
  my $database = #REDACTED#;
  
  
  my $dsn = "DBI:mysql:$database:localhost";
  our $dbh = DBI->connect($dsn,$username,$password);
  }
