#!/usr/bin/perl
# parse_crossmatch_for_Abruijn.pl
# input:  crossmatch result
# output: glue file for RepeatGluer(repeat_sin_read)

# Degui Zhi
# Starts Oct 25, 2004
# Modified Oct 25, 2004

use strict;

# selection criteria
my $SCORE = 350;
my $LENGTH = 40;
my $perc_identity;
my $perc_indel;
my $step = 100;

use lib "/bioinf/oce01/dzhi/bioperl/bioperl-1.0.1"; # on bioinf
use lib '/home/dzhi/bioperl/bioperl-1.0.1';
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long; 
use FileHandle;
use Bio::SearchIO;

my ($blast_result, $name_file);
#get option 
GetOptions ( 
"blast=s"=>\$blast_result,
"Score=s"=>\$SCORE,
"Length=s"=>\$LENGTH,
"name_file=s"=>\$name_file,
);

die "Usage: $0 -b crossmatch_result_file -n seq_names (-L Length_cutoff -S Score_cutoff)  \n" unless ($blast_result and $name_file);


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

print "\nparsing $blast_result\n";

# parse the cross-match output file
open(CM, $blast_result);

while(<CM>)
  {
    last if (/Maximal single base matches \(low complexity regions\):/);
  }

my %turn;
$turn{'query'} = 'hit';
$turn{'hit'} = 'query';

my $turn;
my $hit;
while(<CM>)
  {
    if (/Transitions \/ transversions/)
      {
	if (($hit->{score} >= $SCORE) and ($hit->{length} >= $LENGTH)) {
	  &OutputGlue($hit);
	}
	$num_alignment++;
	next;
      }

    if (/Gap_init rate/)
      {
	next;
      }

    if (/^\s*(\S\s?\S+)\s+(\d+)\s+([A-Z|-]+)\s+(\d+)\s*$/)
      {
	$hit->{seq}{$turn} .= $3;
	$turn = $turn{$turn};
	next;
      }

    if (/^\s*(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\S+).*$/)
      {
	my @items = split /\s+/;
	$hit = ();;
	shift @items unless $items[0]; # to remove spaces at the begining
	$hit->{score} = $items[0];
	$hit->{perc_id} = $items[1];
	$hit->{perc_ins} = $items[2];
	$hit->{perc_del} = $items[3];
	$hit->{query_name} = $items[4];
	$hit->{query_begin} = $items[5];
	$hit->{query_end} = $items[6];
	$hit->{length} = 1+ abs($items[6]-$items[5]);

	if ($items[8] eq 'C')
	  {
	    # need some careful flippings
	    $hit->{strand} = -1;
	    $hit->{hit_name} =  $items[9];
	    $hit->{hit_begin} =  $items[12];
	    $hit->{hit_end} =  $items[11];
	    $hit->{query_begin} =  $items[6];
	    $hit->{query_end} =  $items[5];	    
	  }
	else {
	    $hit->{strand} = 1;
	    $hit->{hit_name} =  $items[8];
	    $hit->{hit_begin} =  $items[9];
	    $hit->{hit_end} =  $items[10];	    
	}

	$turn = 'query';
      }
  }

close HAIXU;
printf "parsed %d alignments\n", $num_alignment;

sub OutputGlue
  {
    my ($hit) = @_;
    
    my $hitindex = $index_of_seq{$hit->{hit_name}};
    my $qindex = $index_of_seq{$hit->{query_name}};
    my $hi = $hit->{hit_begin};
    my $qi = $hit->{query_begin};
#    my $len = $hit->{query_end} - $hit->{query_begin} + 1;
    my $len = length $hit->{seq}{query};
    for my $i (0..$len-1)
      {
	my $sq = substr ($hit->{seq}{query}, $i, 1);
	my $sh = substr ($hit->{seq}{hit}, $i, 1);
	
	if ($sq eq '-')
	  { $hi ++; next; }
	if ($sh eq '-')
	  { $qi += $hit->{strand}; next; }
	
	# $hi+1, $qi+1 just to create artificial unique start for each sequence
	printf HAIXU "%10d %d; \t%10d %d;\n", $hi+1, $hit->{strand}*$hitindex, $qi+1, $qindex; 
	$qi += $hit->{strand};
	$hi ++;
      }
  }
