#!/usr/bin/perl
use lib "."; # put appropriated path to EULER.pm
use EULER;

$dot_file = shift;
$out_file = shift;
$interval_file = shift;
$num_seqs = shift;
$min_edge_length = shift;
$path_file = shift;
$hipath = shift;

# $min_edge_length = 20 unless $min_edge_length;

die "Usage: $0 dot_file out_file interval_file num_seqs min_length [path_file hightlight_path]\n" unless $dot_file and $out_file and $interval_file and $num_seqs;

$out = $out_file;

EULER::LoadEdgeInterval($interval_file, \@intv, \%edge_map);

# insert how to load input_path;
# later do how to translate input path into output path;

my %graph = EULER::LoadDot($dot_file);


foreach $n (keys %{$graph{node}})
  {
    $source{$n} = 1 if (scalar (@{$graph{inedges}{$n}}) == 0);
    $sink{$n}   = 1 if (scalar (@{$graph{outedges}{$n}}) == 0);    
  }

foreach $e (keys %{$graph{edge}})
  {
    $found{$e} = 1 if ($graph{'len'}{$e} > $min_edge_length);
  }

foreach $e (keys %{$graph{edge}})
{
    next if $found{$e};
    $length = &FindSupertangle($e);
}

&WriteSuperGraph;

# END OF MAIN
###############################################################

sub FindSupertangle
{
    my $e = shift;
    $root = $e;
    my $length = 0;

    $found{$e} = 1;

    if (my $from = $graph{edge}{$e}{from}) {  &Traverse($from) ; }
    if (my $to   = $graph{edge}{$e}{to})   {   &Traverse($to); }

    $tangle{$root}{$root} = 1;
    @es =  keys %{$tangle{$e}};
    @{$edges{$root}} = @es;
    foreach $edge (@es)
      {
	$length += $graph{'mul'}{$edge} * $graph{'len'}{$edge};
      }
    print "super-short-edge $e (length=$length):\t", join (",", @es), "\n";
    return $length;
}

sub Traverse
{
    my $node = shift;

    foreach my $e (@{$graph{inedges}{$node}})
    {
	next if ($e == $root);
	next if ($found{$e});

	$found{$e} = 1;
	$tangle{$root}{$e} = 1;
	&Traverse ($graph{edge}{$e}{from});
    }

    foreach my $e (@{$graph{outedges}{$node}})
    {
	next if ($e == $root);
	next if  ($found{$e});

	$found{$e} = 1;
	$tangle{$root}{$e} = 1;
	&Traverse ( $graph{edge}{$e}{to} );
    }
}


sub WriteSuperGraph
  {
open (OUT, ">$out") or die;
print OUT "digraph G {\n";
printf OUT "rankdir=LR;\n";
printf OUT "node [fontsize = 136];\n";
printf OUT "edge [fontsize = 136];\n";
printf OUT "{rank=source; %s}\n", join ";", (keys %source);
printf OUT "{rank=sink; %s}\n", join ";", (keys %sink);
$format = 'pdf';
if ($format eq 'webdot')
{
    print OUT "\tsize=\"18,12\";\n";
}
elsif ($format eq 'pdf')
{
    print OUT "\tsize=\"8,8\";\n";
}

# print node labels
foreach $super (keys %edges)
{
    print OUT "\ts$super [shape=box];\n";
    open (DOT, $dot_file) or die;

    while (<DOT>)
    {
	if (/(\d+) -> (\d+)/)
	{
	    $from = $1;
	    $to = $2;
	}

	if (/label = \"(\d+),(\d+)\((\d+)\)\"/)
	{
	    $edge = $2;
	}

	foreach $e (@{$edges{$super}})
	{
	    if ($e eq $edge)
	    {
		$super{$from} = $super unless $source{$from};
		$super{$to} = $super unless $sink{$to};
	    }
	}
    }
}
foreach $n (keys %{$graph{node}})
  {
    if ($source{$n})
      {
	$e = $graph{outedges}{$n}[0];
	$inp = $edge_map{$e}{intv}[0]{strand}+1;
	$inp = -1*($inp-$num_seqs) if ($inp > $num_seqs);
#	printf OUT "%d [label = \"h%d\", style=\"setlinewidth(10)\"];\n", $n, $inp;
	printf OUT "%d [label = \"%d\", style=\"setlinewidth(10)\"];\n", $n, $inp;
      }
    elsif ($sink{$n})
      {
	$e = $graph{inedges}{$n}[0];
	$inp = $edge_map{$e}{intv}[0]{strand}+1;
	$inp = -1*($inp-$num_seqs) if ($inp > $num_seqs);
	printf OUT "%d [label = \"%d\", style=\"setlinewidth(10)\"];\n", $n, $inp;
      }
    elsif ($super{$n})
      {
	printf OUT "s%d [fontcolor=white, style=\"setlinewidth(10)\"];\n", $super{$n};
      }
    else
      {
	printf OUT "%d [fontcolor=white, style=\"setlinewidth(10)\"];\n", $n;
      }
  }

{
    close DOT;
    open (DOT, $dot_file) or die;
    while (<DOT>)
    {
	if (/(\d+) -> (\d+)/)
	{
	    if ($super{$1} and $super{$2} and ($super{$1} eq $super{$2}) ) {next;}
	    if ($from = $super{$1})
	    {
		s/$1 ->/s$from ->/;
	    }
	    if ($to = $super{$2})
	    {
		s/-> $2/-> s$to/;
	    }
	    s/style=bold,//g;


	if (/(^.*->.*label\s+=\s+)\"(\d+),(\d+)\((\d+)\)\"(.*)/)
	  {
	      $prefix = $1;
	      $len = $2;
	      $edge_id = $3;
	      $mul = $4;
	      $suffix = $5;
	      $style = ", style=\"setlinewidth(10)\", arrowsize = 10";
	      
	      if ($mul > 1)
		{
		  $style = ", style=\"setlinewidth(30)\", arrowsize = 10";
		  $label = "\"$len($mul)\", weight=10,color=red";
		}
	      else
		{
		  $label = "\"$len\"";
		  $label = "\"\"" if $len <= $min_edge_length;
		}
	      $line = $prefix . $label . $style . $suffix . "\n";
	    }
	    
	    print OUT $line;
	    next;
	    
	}
    }
}


print OUT "}\n";

close OUT;

}
