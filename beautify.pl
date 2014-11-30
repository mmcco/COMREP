#!/usr/bin/perl

$intv = shift;
$num_seqs = shift;
$min_len = shift;

die "Usage: $0 INTV NUM_SEQS [MIN_LEN=10] \n" unless $intv;
$min_len = 10 unless $min_len;

for (<*.dot>)
{
  next if /^b\./;
    system("beautify_graph.pl $_ b.$_ $intv $num_seqs $min_len");
#    system("dot -Tps b.$_ -o b.$_.ps");
}
