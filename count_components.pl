#!/usr/bin/perl

# To print a list of connected components of the comparative graph,
# with the breakdown of statistics of each component

# Degui Zhi
# 2004.12 - 2005.1

use lib "."; # put appropriate path to EULER.pm
use EULER;

$super = shift;
$intv = shift;
$ga = shift;
$dir = shift;

die "Usage: $0 SUPERTANGLE INTV GENOME_ASSIGNMENT DIR(=subs)\n" unless $intv;

EULER::LoadEdgeInterval($intv, \@intv, \%edge_map);

my %genome = &LoadGenomeFile($ga);
$num_inputs = ( scalar keys %genome ) / 2;

my @things = qw(id num_edges num_paths num_seqs category num_mouse_compnt num_rat_compnt);
print join " ", @things;
print "\n";

open(SU, $super) or die;
while(<SU>)
  {
    if (/super (\d+):\s*(\S+)/)
      {
	my $super = $1;
	my %in;
	my %compressed_in;
	my %ge_in;
	@th = split ",", $2;
	foreach my $e (@th)
	  {
	    foreach $t (@{$edge_map{$e}{intv}})
	      {
		$in{$$t{strand}} = 1;
		$compressed_in{$$t{strand}%$num_inputs} = 1;
	      }
	  }

	foreach (keys %in)
	  {
	    $ge_in{$genome{$_}} ++;
	  }
	my @in_ges = keys %ge_in;
	my $cat;
	if ( (scalar @in_ges) == 1)
	  {
	    $cat = $in_ges[0];
	  }
	else
	  {
	    $cat = 2;
	  }

	my $edge_count = 0;
	$dir ||= "subs";
	my $subfile = "$dir/subgraph_s$super.dot";
	if (-e $subfile) 
	  {
	    open(SF, $subfile) or die;
	    while(<SF>)
	      {
		$edge_count ++ if /->/;
	      }
	    close SF;
	  }
	else
	  { $edge_count = 1; }

	printf "%d %d %d %d %d %d %d\n", $super, $edge_count, scalar keys %compressed_in, scalar keys %in, $cat, $ge_in{0}, $ge_in{1};
      }
  }

sub LoadGenomeFile
  {
    my $file = shift;
    my %genome;

    if (not -e $file)
      {
	for my $i (0..10000)
	  {
	    $genome{$i} = 1;
	  }
	return %genome;
      }


    open(FI, $file) or die;
    while(<FI>)
      {
	next if /^\#/;
	chomp;
	@line = split /\s+/;
	for my $i ($line[1]..$line[2])
	  {
	    $genome{$i-1} = $line[0];
	  }
      }
    close FI;
    return %genome;
  }
