#!/usr/bin/perl
use lib "."; # put appropriate path to EULER.pm
use EULER;

$seqname = shift;
$interval_file = shift;
$min_edge_length = shift;


die "Usage: $0 seqname-file interval_file min_length\n" unless $seqname and  $interval_file and $min_edge_length;

EULER::LoadEdgeInterval($interval_file, \@intv, \%edge_map);
my @name;
my %seq_length;
$num_seqs=&LoadSeqNames($seqname, \@name, \%seq_length); # names are duplicated

foreach $e (keys %edge_map)
  {
    next if $edge_map{$e}{length} <= $min_edge_length;
    $num_intv = scalar @{$edge_map{$e}{intv}};
    my $done = 0;
    for ($i=0; $i<$num_intv; $i++)
      {
	$n1 = $name[${$edge_map{$e}{intv}}[$i]{strand}];
	for ($j = 0; $j < $i; $j++)
	  {
	    $n2 = $name[${$edge_map{$e}{intv}}[$j]{strand}];
	    if (&dist($n1,$n2))
	      {
#		printf "edge %d(%d,%d): $n1, $n2\n", $e, 
#		   $edge_map{$e}{length},  $edge_map{$e}{mul};
		$done = 1;
		last;
	      }
	  }
	last if $done;
      }

    # collect stats of repeats passing through this edge
    if ($done)
      {
	for ($i=0; $i<$num_intv; $i++)
	  {
	    next if ${$edge_map{$e}{intv}}[$i]{strand} >= $num_seqs;
	    $n1 = $name[${$edge_map{$e}{intv}}[$i]{strand}];
	    $stats{$n1}{num_good_edges} ++;
	    $stats{$n1}{len_good_edges} += $edge_map{$e}{length};
	  }
      }
  }

@ids = reverse sort {$stats{$a}{num_good_edges}<=> $stats{$b}{num_good_edges}} keys %stats;

foreach $sid (@ids)
  {
    printf "%20s %10d %10d %8.2f\n", 
      $sid, $stats{$sid}{num_good_edges},
	$stats{$sid}{len_good_edges},
	  $stats{$sid}{len_good_edges}/$seq_length{$sid};
  }

sub LoadSeqNames
  {
    my ($namefile, $name_, $seq_length_) = @_;

    my $num_seqs;
    @$name_ = ();
    open(NA, $namefile) or die;
    while(<NA>)
      {
	if (/>(\S+)/)
	  {
	    $name = $1;
	    push @$name_, $name;
	  }
	else
	  {
	    chomp;
	    $$seq_length_{$name} += length($_);
	  }
      }
    close NA;

    # duplicate it
    $num_seqs = scalar @$name_;
    for ($i=0; $i<$num_seqs; $i++)
      {
	push @$name_, $$name_[$i];
      }
    return $num_seqs;
  }

sub dist
  {
    my ($n1, $n2) = @_;
    if (substr($n1,0,2) eq substr($n2,0,2)) {return 0;}
    else {return 1};
  }



