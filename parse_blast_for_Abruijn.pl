#!/usr/bin/perl
# parse_blast_for_Abruijn.pl
# input:  blastall -p blastp result
# output: glue file for RepeatGluer(repeat_sin_read)

# Degui Zhi
# Starts Feb 03, 2004
# Modified Sep 29, 2004

use strict;

my $LENGTH = 20;
my $FRAC = 0.1; # frac_conserved
my $E_value = -1;
my $step = 100;


use lib "/bioinf/oce01/dzhi/bioperl/bioperl-1.0.1"; # on bioinf
use lib '/home/dzhi/bioperl/bioperl-1.0.1';
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long; 
use FileHandle;
use Bio::SearchIO;

my ($blast_result_list, $blast_result, $name_file);
#get option 
GetOptions ( 
"blast=s"=>\$blast_result,
"Frac=s"=>\$FRAC,
"E-value=s"=>\$E_value,
"Length=s"=>\$LENGTH,
"name_file=s"=>\$name_file,
);

die "Usage: $0 -b blast_result_file -n seq_names (-F min_Fraction_conserved) (-L MIN_LENGTH) (-E min_E-value) \n" unless ($blast_result and $name_file);


my $num_alignment = 0;

my @blast_results;
@blast_results = ($blast_result);

my @seq_names;
my %index_of_seq;
open (NAME, $name_file) or die;
while (<NAME>)
  {
    if (/>(\S+)\s/)
      {
	push @seq_names, $1;
	$index_of_seq{$1} = (scalar @seq_names) ;
      }
  }
close NAME;

my $haixu = $blast_result.".map"; # only record mappings for repeat reads

open (HAIXU, ">$haixu") or die;
print HAIXU "1000000\n"; # as required by repeat_sin_read
my $count = 0;

    print "\nparsing $blast_result\n";
    my $searchio = new Bio::SearchIO (-format => 'blast',
				      -file   => $blast_result);

while (my $blast = $searchio->next_result)
  {
    my $qid = $blast->query_name;
    my $qindex = $index_of_seq{$qid};
    next unless $qindex;
    while (my $hit = $blast->next_hit)
      {
#	my $hitid = $hit->hit_name;
	my $hitid = $hit->name;
	if ($qid eq $hitid)
	  {
	    print "je!\n";
	  }
	my $hitindex = $index_of_seq{$hitid};
	next unless $hitindex;
	while( my $hsp = $hit->next_hsp )
	  {
	    if ($E_value > 0)
	      {
		next unless ( ($hsp->evalue() < $E_value));
	      }
	    else
	      {
		next unless ( ($hsp->length('total') > $LENGTH) and (($hsp->num_conserved > $FRAC * $hsp->length('total')) or ( $hsp->num_identical > $FRAC * $hsp->length('total') )) );
	      }
	    $num_alignment ++;
	    if ($hsp->strand('hit') >= 0) {
	      my $hi = $hsp->start('hit');
	      my $qi = $hsp->start('query');
	      for my $i (0..$hsp->length('total')-1)
		{
		  my $sq = substr ($hsp->query_string, $i, 1);
		  my $sh = substr ($hsp->hit_string, $i, 1);
		  
		  if ($sq eq '-')
		    { $hi ++; next; }
		  if ($sh eq '-')
		    { $qi ++; next; }
		  
		  # $hi+1, $qi+1 just to create artificial unique start for each sequence
		  printf HAIXU "%10d %d; \t%10d %d;\n", $hi+1, $hitindex, $qi+1, $qindex; 
		  $hi ++;
		  $qi ++;
		}
	    } else {
	      my $hi = $hsp->end('hit');
	      my $qi = $hsp->start('query');
	      for my $i (0..$hsp->length('total')-1)
		{
		  my $sq = substr ($hsp->query_string, $i, 1);
		  my $sh = substr ($hsp->hit_string, $i, 1);
		  
		  if ($sq eq '-')
		    { $hi --; next; }
		  if ($sh eq '-')
		    { $qi ++; next; }
		  
		  # $hi+1, $qi+1 just to create artificial unique start for each sequence
		  printf HAIXU "%10d -%d; \t%10d %d;\n", $hi+1, $hitindex, $qi+1, $qindex; 
		  $hi --;
		  $qi ++;
		}
	    }
	  }
      }
  }

close HAIXU;
printf "parsed %d alignments\n", $num_alignment;

