#!/usr/bin/perl

use Getopt::Long;


my $program_name = $0;
GetOptions (
	     "subs=s"=>\$subs_file,
	     "dot=s"=>\$dot_file,
	     "format=s"=>\$format,
	     "out=s"=>\$out,
	   );

die "Usage: $program_name -s subgraph_file -d dot_file [-f output_format=webdot] [-o out_file]\n" unless ( ($subs_file and $dot_file) );
$format='webdot' unless $format;
$out = "super.dot" unless $out;

open (SUBS, $subs_file) or die;
while (<SUBS>)
{
    if (/super (\d+):\s+(.*)/)
    {
	@{$edges{$1}} = split ",", $2;
    }
}


open (OUT, ">$out") or die;
print OUT "digraph G {\n";

if ($format eq 'webdot')
{
    print OUT "\tsize=\"18,12\";\n";
    
    print OUT "node [URL=\"http://bluegene.ucsd.edu/cgi-bin/webdot/nm_dot/subgraph_\\N.dot.gif\"];\n";
}
elsif ($format eq 'pdf')
{
    print OUT "\tsize=\"8,8\";\n";
    print OUT "node [URL=\"subgraph_\\N.dot.gif\"];\n";   
}

foreach $super (keys %edges)
{

    print OUT "\ts$super [shape=box,color=red];\n";
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
		$super{$from} = $super ;
		$super{$to} = $super ;
	    }
	}
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
	    print OUT;
	    next;
	    
	}
    }
}


print OUT "}\n";
close OUT;
