#!/usr/bin/perl -w
##Eric Smith

use strict;
use DBI;
prepareDbh();

my $coverage_table = 'genomes_TAPL_coverage';

$dbh->do(qq{drop table if exists $coverage_table});
$dbh->do(qq{create table $coverage_table
            (gtc_id int auto_increment primary key,
			 sample varchar(20),
			 species varchar(10),
			 sex varchar(1),
			 country varchar(20),
			 region varchar(20),
			 pos int,
			 cov int,
			 index(sample))});

my %sample_lookup;
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
		}
	  $j++;
	  }
	}
  else
    {
	my @values = split /\t/, $line;
	my $position = $values[1];
	if ($position < 11400000 and $position > 11100000)
	  {
	  my $k = 9;
	  while ($k < scalar @values)
	    {
	    my $sample_id = $sample_lookup[$k];
	    my $sample_allele;
		my @gt_field = split /:/, $values[$k];
		my $genotype = $gt_field[0];
	    unless ($genotype eq './.')
		  {
		  my $read_depth = $gt_field[2];
		  if ($read_depth =~ /\,/)
		    {
		    $read_depth = $gt_field[3];	
		    }
		  my ($form, $sex, $country, $region) = $dbh->selectrow_array(qq{select form, sex, country, region
		                                                                 from lk_genomes_project_samples
																		 where sample_id like '${sample_id}%'});
		  
		  my $species;
		  if ($form eq 'M')
		    {
		    $species = 'coluzzii';	
		    }
		  else
		    {
		    $species = 'gambiae';	
		    }
		  $dbh->do(qq{insert into $coverage_table
			          values
					  (null, '$sample_id', '$species', '$sex', '$country', '$region', $position, $read_depth)});
		  }
		$k++;
		}
	  }	
	} 	
  }

sub prepareDbh
  {
  my $username = "david";
  my $password = "drosophila";
  my $database = "tapl";
  
  
  my $dsn = "DBI:mysql:$database:localhost";
  our $dbh = DBI->connect($dsn,$username,$password);
  }
