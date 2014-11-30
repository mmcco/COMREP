#!/usr/bin/perl

use Getopt::Long;


my $program_name = $0;
GetOptions (
	     "subs=s"=>\$subs_file,
	     "dot=s"=>\$dot_file,
	   );

die "Usage: $program_name -s subgraph_file -d dot_file \n" unless ( ($subs_file and $dot_file) );


open (SUBS, $subs_file) or die;
while (<SUBS>)
{
    if (/super (\d+):\s+(.*)/)
    {
	@{$edges{$1}} = split ",", $2;
    }
}

foreach $super (keys %edges)
{

    open (OUT, ">subgraph_s$super.dot") or die;
    print OUT "digraph G {\n";
    print OUT "\tsize=\"8,8\";\n";
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
		$nodes{$super}{$from} = 1 ;
		$nodes{$super}{$to} = 1 ;
	    }
	}
    }

    close DOT;
    open (DOT, $dot_file) or die;
    while (<DOT>)
    {
	if (/(\d+) -> (\d+)/)
	{
	    print OUT if ($nodes{$super}{$1} or $nodes{$super}{$2});
	}
    }
   

    print OUT "}\n";

    close OUT;
}
