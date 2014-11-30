#!/usr/bin/perl
$dot = shift;
$min_len = shift;
$min_mul = shift;
die "Usage: $0 DOT(beautified) min_len min_mul \n" unless $dot and $min_mul;

$eid = 1;

open(DOT, $dot);
while(<DOT>)
  {
    if ((/(\S+) -> (\S+) \[label = \"(.+)\"/))
      {
	$from = $1;
	$to = $2;
	$label = $3;
	if ($label =~/(\d+)\((\d+)\)/)
	  {
	    $len = $1;
	    $mul = 2;
	  }
	elsif ($label =~/(\d+)/)
	  {
	    $len = $1;
	    $mul = 1;
	  }
	else
	  {
	    $len = 0;
	    $mul = 1;
	  }
	next if $len<$min_len;
	next if $mul<$min_mul;
	
	$rg{edge}{$eid}{len} = $len;
	$rg{edge}{$eid}{mul} = $mul;
	$rg{edge}{$eid}{from} = $from;
	$rg{edge}{$eid}{to} = $to;
	$rg{node}{$from}{outedges}{$eid} = 1;
	$rg{node}{$to}{inedges}{$eid} = 1;
	$eid ++;
	next;
      }
    # for un-beautified graphs
    if ((/(\S+) -> (\S+) \[label = \"(\d+),\d+\((\d+)\)\"/))
      {
	next if $3<$min_len or $4<$min_mul;
	
	$rg{edge}{$eid}{len} = $3;
	$rg{edge}{$eid}{mul} = $4;
	$rg{edge}{$eid}{from} = $1;
	$rg{edge}{$eid}{to} = $2;
	$rg{node}{$1}{outedges}{$eid} = 1;
	$rg{node}{$2}{inedges}{$eid} = 1;
	$eid ++;
      }
  }

# finding knot nodes
foreach $n (keys %{$rg{node}})
  {
    my $ins = scalar keys %{$rg{node}{$n}{inedges}};
    my $outs = scalar keys %{$rg{node}{$n}{outedges}};
    if (($ins > 0) and ($outs > 0) and ($ins+$outs>2))
      {
	$knot{$n} = 1;
      }
  }

#output knots
return unless scalar keys %knot;
print "digraph G{\n";
foreach $n (keys %knot)
  {
    print "\t$n [shape=box]\n";
    foreach $e (keys %{$rg{node}{$n}{inedges}})
      {
	printf "\t$n%s -> $n%s [label=\"%d(%d)\"];\n",
	  $rg{edge}{$e}{from}, $rg{edge}{$e}{to},
	    $rg{edge}{$e}{len}, $rg{edge}{$e}{mul};
      }
    foreach $e (keys %{$rg{node}{$n}{outedges}})
      {
	printf "\t$n%s -> $n%s [label=\"%d(%d)\"];\n",
	  $rg{edge}{$e}{from}, $rg{edge}{$e}{to},
	    $rg{edge}{$e}{len}, $rg{edge}{$e}{mul};
      }
  }
print "}\n";
