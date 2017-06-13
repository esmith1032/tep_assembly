#!/usr/bin/perl -w
##Eric Smith

use strict;
use DBI;
prepareDbh();

$dbh->do(qq{drop table if exists haplotype_consensus_alns_snps});
$dbh->do(qq{create table haplotype_consensus_alns_snps
            (hcas_id int auto_increment primary key,
	           pest_pos int,
			       align_pos int,
             pest_nt varchar(1),
             Ra_nt varchar(1),
             Rb_nt varchar(1),
             S_nt varchar(1),
             index(pest_pos))});

$dbh->do(qq{drop table if exists haplotype_consensus_alns_indels});
$dbh->do(qq{create table haplotype_consensus_alns_indels
            (hcai_id int auto_increment primary key,
	           pest_pos int,
			       align_pos int,
             len int,
             pest_has varchar(10),
			       Ra_has varchar(10),
			       Rb_has varchar(10),
			       S_has varchar(10),
             sequence longtext,
             index(pest_pos))});
			 
$dbh->do(qq{drop table if exists haplotype_consensus_alns_seqs});
$dbh->do(qq{create table haplotype_consensus_alns_seqs
			(haplotype varchar(10),
			 seq longtext,
			 index(haplotype))});
			 
$dbh->do(qq{drop table if exists haplotype_consensus_alns_ns});
$dbh->do(qq{create table haplotype_consensus_alns_ns
            (haplotype varchar(5),
			       raw_pos int,
			       align_pos int,
			       n_len int,
			       index(haplotype))});

my $pest_start = 11079073;
my $align_file = '/Research/gambiae/tapl/spades/all_bacs/pacbio_hybrid/Ra_Rb_PEST_Aligns.fasta';

my %seq_info;

open ALIGNS, "$align_file";
while (my $line = <ALIGNS>) {
  chomp $line;
  if ($line =~ /^>/) {
    my $seq_name = substr($line, 1, length($line) - 1);
	  my $seq = <ALIGNS>;
	  chomp $seq; 
    $seq_info{$seq_name}{seq} = $seq;
	  $seq_info{$seq_name}{len} = length($seq);
	  $seq_info{$seq_name}{indel} = 'N';
	  $seq =~ m/[ACTGN]/g;
	
	  $seq_info{$seq_name}{start} = pos($seq) - 1;
    $dbh->do(qq{insert into haplotype_consensus_alns_seqs
				values
				('$seq_name', '$seq')});
	  my $k = 0;
    my $raw_pos = 0;
	  while ($k <= length($seq)) {
      my $nt = substr($seq, $k, 1);
	    if ($nt eq 'N') {
	      my $new_seq = substr($seq, $k, length($seq) - $k);
		    $new_seq =~ m/[ACTG\-]/g;
		    my $n_len = pos($new_seq);
		    $dbh->do(qq{insert into haplotype_consensus_alns_ns
		            values
					      ('$seq_name', $raw_pos, $k, $n_len)});
		    $k += $n_len;
		    $raw_pos += $n_len;
		  }
	    
      unless ($nt eq '-') {
	  	  $raw_pos++;
	    }
	    
      $k++;
	  }
	}
}
close ALIGNS;  

my $len;

LOOP: foreach my $key (sort {$seq_info{$b}{len} <=> $seq_info{$a}{len}} keys %seq_info) {
  $len = $seq_info{$key}{len};
  last LOOP;
}

my $pest_offset = 0;
my $i = 0;
SNP_LOOP: while ($i <= $len) {
  my %nts;
  foreach my $haplotype (keys %seq_info) {
    if ($i <= $seq_info{$haplotype}{len} and $i >= $seq_info{$haplotype}{start}) {
      my $base = substr($seq_info{$haplotype}{seq}, $i, 1);
	    $nts{$haplotype} = $base;
    } else {
	    $nts{$haplotype} = 'x';	
	  }	  
  }
  my $hap_1 = 'PEST';	

  if ($nts{$hap_1} eq '-' or $nts{$hap_1} eq 'N') {
	  foreach my $hap_2 (keys %nts) {
	    next if ($hap_2 eq $hap_1);
      next if ($nts{$hap_2} eq '-');
	    next if ($nts{$hap_2} eq 'N');
	    next if ($nts{$hap_2} eq 'x');
	  
	    foreach my $hap_3 (keys %nts) {
	      next if ($hap_3 eq $hap_2);
		    next if ($hap_3 eq $hap_1);
		    next if ($nts{$hap_3} eq '-');
		    next if ($nts{$hap_3} eq 'N');
		    next if ($nts{$hap_3} eq 'x');
		
		    if ($nts{$hap_2} ne $nts{$hap_3}) {
	        $dbh->do(qq{insert into haplotype_consensus_alns_snps
		                  values
					           (null, $pest_start + $i - $pest_offset, $i, if('$nts{PEST}' = 'x', null, '$nts{PEST}'), if('$nts{Ra}' = 'x', null, '$nts{Ra}'), if('$nts{Rb}' = 'x', null, '$nts{Rb}'), if('$nts{S}' = 'x', null, '$nts{S}'))});
		  
		      if ($nts{$hap_1} eq '-') {
            $pest_offset += 1;
          }
		      $i++;
		      next SNP_LOOP;	
		    }
	    }
	  }
  } else {
	  foreach my $hap_2 (keys %nts) {
      next if ($hap_1 eq $hap_2);
      next if ($nts{$hap_2} eq 'x');
      next if ($nts{$hap_2} eq 'N');
      next if ($nts{$hap_2} eq '-');		
    
      if ($nts{$hap_1} ne $nts{$hap_2}) {
    	  $dbh->do(qq{insert into haplotype_consensus_alns_snps
    	              values
    	              (null, $pest_start + $i - $pest_offset, $i, if('$nts{PEST}' = 'x', null, '$nts{PEST}'), if('$nts{Ra}' = 'x', null, '$nts{Ra}'), if('$nts{Rb}' = 'x', null, '$nts{Rb}'), if('$nts{S}' = 'x', null, '$nts{S}'))});	
		    $i++;	
		    next SNP_LOOP;
		  }	
    }
  }		
  $i++;
} 

$pest_offset = 0;
$i = 0;
INDEL_LOOP: while ($i <= $len) {
  my $pest_nt;
  my $Ra_nt;
  my $Rb_nt;
  my $S_nt;
  
  if ($i >= $seq_info{PEST}{start} and $i < $seq_info{PEST}{len}) {
    $pest_nt = substr($seq_info{PEST}{seq}, $i, 1);
	} elsif ($i < $seq_info{PEST}{start}) {
    $pest_nt = 'S';
  } else {
    $pest_nt = 'E';
	}
	
  if ($i >= $seq_info{Ra}{start} and $i < $seq_info{Ra}{len}) {
    $Ra_nt = substr($seq_info{Ra}{seq}, $i, 1);
  } elsif ($i < $seq_info{Ra}{start}) {
    $Ra_nt = 'S';
  } else {
    $Ra_nt = 'E';
	}	
	
  if ($i >= $seq_info{Rb}{start} and $i < $seq_info{Rb}{len}) {
    $Rb_nt = substr($seq_info{Rb}{seq}, $i, 1);
  } elsif ($i < $seq_info{Rb}{start}) {
    $Rb_nt = 'S';
  }	else {
    $Rb_nt = 'E';
  }

  if ($i >= $seq_info{S}{start} and $i < $seq_info{S}{len}) {
    $S_nt = substr($seq_info{S}{seq}, $i, 1);
  } elsif ($i < $seq_info{S}{start}) {
    $S_nt = 'S';
	}	else {
    $S_nt = 'E';
  }
  
  if ($pest_nt eq '-' or $Ra_nt eq '-' or $Rb_nt eq '-' or $S_nt eq '-') {
    my $pest_has;
	  my $Ra_has;
	  my $Rb_has;
	  my $S_has;
	  my $pest_seq_len = $len;
	  my $Ra_seq_len = $len;
	  my $Rb_seq_len = $len;
	  my $S_seq_len = $len;
	  my $indel_len = $len;
	  my $found_seq = 'N';
	  my $indel_seq;
	  my $n_len = $len;
	  
	  if ($pest_nt eq '-') {
	    $pest_has = 'N';
      my $new_nt_seq = substr($seq_info{PEST}{seq}, $i, $seq_info{PEST}{len} - $i);
	    $new_nt_seq =~ m/[ACTGN]/g;
	    my $new_seq_start = pos($new_nt_seq);
	    if ($new_seq_start < $indel_len) {
        $indel_len = $new_seq_start - 1;
      }
	  } elsif ($pest_nt eq 'A' or $pest_nt eq 'C' or $pest_nt eq 'T' or $pest_nt eq 'G' or $pest_nt eq 'N') {
	    $pest_has = 'Y';	
      my $new_indel_seq = substr($seq_info{PEST}{seq}, $i, $seq_info{PEST}{len} - $i);
	  
      if ($new_indel_seq =~ m/\-/g) {
        $pest_seq_len = pos($new_indel_seq);
      } else {
        $pest_seq_len = $seq_info{PEST}{len} - $i;
      }	
	  
      if ($pest_seq_len < $indel_len) {
	      $indel_len = $pest_seq_len - 1;
		  }
	  } elsif ($pest_nt eq 'S') {
      $pest_has = 'start';
	  } elsif ($pest_nt eq 'E') {
      $pest_has = 'end';
    }      

  	if ($Ra_nt eq '-') {
  	  $Ra_has = 'N';
      my $new_nt_seq = substr($seq_info{Ra}{seq}, $i, $seq_info{Ra}{len} - $i);
  	  $new_nt_seq =~ m/[ACTGN]/g;
  	  my $new_seq_start = pos($new_nt_seq);
  	  
      if ($new_seq_start < $indel_len) {
        $indel_len = $new_seq_start - 1;
      }
  	} elsif ($Ra_nt eq 'A' or $Ra_nt eq 'C' or $Ra_nt eq 'T' or $Ra_nt eq 'G' or $Ra_nt eq 'N') {
  	  $Ra_has = 'Y';	
      my $new_indel_seq = substr($seq_info{Ra}{seq}, $i, $seq_info{Ra}{len} - $i);
  	  $new_indel_seq =~ m/\-/g;
	    
      if ($new_indel_seq =~ m/\-/g) {
        $Ra_seq_len = pos($new_indel_seq);
      } else {
        $Ra_seq_len = $seq_info{Ra}{len} - $i;
      }	
  	  
      if ($Ra_seq_len < $indel_len) {
  	    $indel_len = $Ra_seq_len - 1;
  		}
  	} elsif ($Ra_nt eq 'S') {
      $Ra_has = 'start';
  	} elsif ($Ra_nt eq 'E') {
      $Ra_has = 'end';
    }      

  	if ($Rb_nt eq '-') {
  	  $Rb_has = 'N';
      my $new_nt_seq = substr($seq_info{Rb}{seq}, $i, $seq_info{Rb}{len} - $i);
  	  $new_nt_seq =~ m/[ACTGN]/g;
  	  my $new_seq_start = pos($new_nt_seq);
  	  
      if ($new_seq_start < $indel_len) {
        $indel_len = $new_seq_start - 1;
      }
  	} elsif ($Rb_nt eq 'A' or $Rb_nt eq 'C' or $Rb_nt eq 'T' or $Rb_nt eq 'G' or $Rb_nt eq 'N') {
  	  $Rb_has = 'Y';	
      my $new_indel_seq = substr($seq_info{Rb}{seq}, $i, $seq_info{Rb}{len} - $i);
  	  $new_indel_seq =~ m/\-/g;
	  
      if ($new_indel_seq =~ m/\-/g) {
        $Rb_seq_len = pos($new_indel_seq);
      } else {
        $Rb_seq_len = $seq_info{Rb}{len} - $i;
      }	
  	  
      if ($Rb_seq_len < $indel_len) {
  	    $indel_len = $Rb_seq_len - 1;
  		}
  	} elsif ($Rb_nt eq 'S') {
      $Rb_has = 'start';
  	} elsif ($Rb_nt eq 'E') {
      $Rb_has = 'end';
    }      

  	if ($S_nt eq '-') {
  	  $S_has = 'N';
      my $new_nt_seq = substr($seq_info{S}{seq}, $i, $seq_info{S}{len} - $i);
  	  $new_nt_seq =~ m/[ACTGN]/g;
  	  my $new_seq_start = pos($new_nt_seq);

  	  if ($new_seq_start < $indel_len) {
        $indel_len = $new_seq_start - 1;
      }
  	} elsif ($S_nt eq 'A' or $S_nt eq 'C' or $S_nt eq 'T' or $S_nt eq 'G' or $S_nt eq 'N') {
  	  $S_has = 'Y';	
      my $new_indel_seq = substr($seq_info{S}{seq}, $i, $seq_info{S}{len} - $i);
  	  $new_indel_seq =~ m/\-/g;
	  
      if ($new_indel_seq =~ m/\-/g) {
        $S_seq_len = pos($new_indel_seq);
      } else {
        $S_seq_len = $seq_info{S}{len} - $i;
      }	
  	  
      if ($S_seq_len < $indel_len) {
  	    $indel_len = $S_seq_len - 1;
  		}
  	} elsif ($S_nt eq 'S') {
      $S_has = 'start';
  	} elsif ($S_nt eq 'E') {
      $S_has = 'end';
    }      
    
	  if ($pest_has eq 'Y' and $pest_nt ne 'N') {
	    my $n_num = () = substr($seq_info{PEST}{seq}, $i, $indel_len) =~ /N/;
      if ($n_num < $n_len) {
  	    $indel_seq = substr($seq_info{PEST}{seq}, $i, $indel_len);
		    $n_len = $n_num;	
	    }
	  }
	
    if ($Ra_has eq 'Y' and $Ra_nt ne 'N') {
      my $n_num = () = substr($seq_info{Ra}{seq}, $i, $indel_len) =~ /N/;
  	  if ($n_num < $n_len) {
    	  $indel_seq = substr($seq_info{Ra}{seq}, $i, $indel_len);
  		  $n_len = $n_num;	
  	  }
    }
	
    if ($Rb_has eq 'Y' and $Rb_nt ne 'N') {
      my $n_num = () = substr($seq_info{Rb}{seq}, $i, $indel_len) =~ /N/;
      if ($n_num < $n_len) {
    	  $indel_seq = substr($seq_info{Rb}{seq}, $i, $indel_len);
  		  $n_len = $n_num;	
  	  }
    }
	
    if ($S_has eq 'Y' and $S_nt ne 'N') {
      my $n_num = () = substr($seq_info{S}{seq}, $i, $indel_len) =~ /N/;
  	  if ($n_num < $n_len) {
        $indel_seq = substr($seq_info{S}{seq}, $i, $indel_len);
  		  $n_len = $n_num;	
  	  }
    }

	  $dbh->do(qq{insert into haplotype_consensus_alns_indels
	              values
			  	      (null, $pest_start + $i - $pest_offset, $i, $indel_len, '$pest_has', '$Ra_has', '$Rb_has', '$S_has', '$indel_seq')});
	
	  if ($pest_has eq 'N') {
      $pest_offset += $indel_len;
	  }
	
    $i += $indel_len;
	  next INDEL_LOOP;
  }	
  $i++;
}

print STDOUT "Finished\n";  

sub prepareDbh {
  my $username = #REDACTED#;
  my $password = #REDACTED#;
  my $database = #REDACTED#;
  
  my $dsn = "DBI:mysql:$database:localhost";
  our $dbh = DBI->connect($dsn,$username,$password);
}
