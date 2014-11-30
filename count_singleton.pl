#!/usr/bin/perl
use lib "."; # put appropriate path to EULER.pm
use EULER;

$file = shift;
%g = EULER::LoadDot($file) or die;

foreach $n (keys %{$g{node}}) 
  {
    $source{$n} = 1 if (scalar (@{$g{inedges}{$n}}) == 0);
    $sink{$n}   = 1 if (scalar (@{$g{outedges}{$n}}) == 0);    
  }


foreach $e (keys %{$g{edge}}) 
  {
    $num_singletons ++ if (($source{ $g{edge}{$e}{from} }) and ($sink{ $g{edge}{$e}{to} }));
  }

print $num_singletons, "\n";

