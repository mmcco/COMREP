######################################################################
# EULER.pm   package for parsing EULER outputs
# Version 1
# Author: Degui Zhi
# Date:   Sep 16, 2004
# Description: Merge from TR_LOAD and NEW_LOAD
#
############################################################################
# use lib "/usr/local/bioperl-1.0";
package EULER;

use 5.005;
use vars qw(@ISA);
use strict;
# use warnings;
use lib "/home/dzhi/bioperl/bioperl-1.0.1";
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::Blast;
use Bio::Tools::Blast::HSP;
use Getopt::Long; 
# use Carp::Assert;


require Exporter;

@ISA = qw(Exporter);


################################################
# sub LoadPathFile
################################################

sub LoadPathFile
  {
    my ($path_file) = @_;
    open(PATH, $path_file) or die;

    my %path;
    my @path;
    my $sid;

    while (<PATH>)
      {
	my $line = $_;
	if ((/^SEQ(\d+) (.*)/) or (/Sequence(\d+)/))
	  {
	    $sid = $1;
	    @{$path{$sid}{chain}} = ();
	  }

	@path = split /\s+/, $line;
	
	while ((my $ele = shift @path) ne undef)
	  {
	    if ($ele =~ /(\d+)\((\d+),(\d+)\)/)
	      {
		my $e = $1+1;
		$path{$sid}{edge}{$e}{mul} = $2;
		$path{$sid}{edge}{$e}{len} = $3;
		$path{$sid}{index}{$e} = scalar @{$path{$sid}{chain}};
  		push @{$path{$sid}{chain}}, $e;
	      }
	  }
	
      }

    return %path;
  }

################################################
# sub LoadMap
################################################

sub LoadMap
{
    my ($mapping_file,$edges_ref ) = @_;
    my %isedge;
    my %map;
    my %seen;
    foreach my $e (@$edges_ref)
    {
	$isedge{$e} = 1;
    }

    open (MAP, $mapping_file) or die;
    while (<MAP>)
    {
	my ($index, $id, $num, $a, $z, $pa, $pz, $path) = split;
	# next unless ($isedge{$a} or $isedge{$z});
	print "warning dup id: $id\n" if $seen{$id};
	$seen{$id} ++;
	$map{$id}{index} = $index;
	$map{$id}{start_edge} = $a;
	$map{$id}{end_edge} = $z;
	$map{$id}{start_pos} = $pa;
	$map{$id}{end_pos} = $pz;
	@{$map{$id}{path}} = split ",", $path;
    }
    return %map;
    
}



################################################
# sub DumpAllEdges
################################################
sub DumpAllEdges
  {
      my $edge_hash_ref = shift;
      foreach my $edge (keys %{$edge_hash_ref})
      {
	  my $outfile = 'edge_' . $edge . '.fa';
	  my $out  = Bio::SeqIO->new('-file' => ">$outfile",
				     '-format' => 'Fasta');
	  $out->write_seq($$edge_hash_ref{$edge});
      }
  }



################################################
# sub LoadEdges
################################################
sub LoadEdges
  {
    my ($edge_file,$verbose,$debug ) = @_;
    my $edgeseq;
    my %edge_hash;

    my $edge_io  = Bio::SeqIO->new('-file' => "$edge_file",
				   '-format' => 'Fasta');

    while ( $edgeseq = $edge_io->next_seq() ) 
    {
	    $edge_hash{$edgeseq->id} = $edgeseq;
    }
    print "[TR_LOAD], # edge loaded: ", scalar (keys %edge_hash), "\n" if $verbose;
   return %edge_hash;
  }



################################################
# sub LoadGraph
################################################

=item &LoadGraph
Input:  $graph_file
Output: %graph_obj
=cut

sub LoadGraph
  {

    my ($graph_file, $verbose, $debug, $option) = @_;
     my $i;
    my ($watson, $crick); 
    my $vertex_no;
    my $no_of_out;
    my $no_of_in;
    my $seq_length;

    my %graph_obj;

    print "in LoadGraph: option = $option\n" if $verbose;
    open (GRAPH, $graph_file) or die "can not open graph file\n";
    
    while (<GRAPH>)
      {
	if (/^Vertex (\d+) (\d+) (\d+) (\d+)/)
	  {
	    $vertex_no = $1+1;
	    $vertex_no = $1 if $option eq 2;
	    $no_of_out = $2;
	    $no_of_in = $3;
	    $seq_length = $4;
	    $graph_obj{node}{$vertex_no} = 1;
	    next;
	  }
	elsif (/^Vertex (\d+) (\d+) (\d+)/)
	  {
	    $vertex_no = $1+1;
	    $vertex_no = $1 if $option eq 2;
	    $no_of_out = $2;
	    $no_of_in = $3;
	    $graph_obj{node}{$vertex_no} = 1;
	    next;
	  }
	
	if (/^SEQ (.+)/){next;}
	  
	if (/Last_edge (.*)$/)
	  {
	    my @aa = split /\s+/, $1;
	    if ($option eq 1) {
		for ($i = 0; $i < scalar @aa; $i ++) { $aa[$i] ++;} # adjust to get the edge number
	    }
	    my %bb = @aa;
	    @{$graph_obj{inedges}{$vertex_no}} = keys %bb;
	    
	    for ($i = 0; $i < scalar @aa; $i += 2)
	      {
		# do isomorphism assignment of edges;
		$watson = $aa[$i];
		$crick  = $aa[$i+1];
#		$watson = $aa[$i] +1 if $option eq 2;
#		$crick  = $aa[$i+1]+1 if $option eq 2;
		$graph_obj{pal}{$watson} = $crick;
		$graph_obj{pal}{$crick} = $watson;
		
		$graph_obj{edge}{$watson}{to} = $vertex_no;
	      }
	  }
 	
	if (/Next_edge (.*)$/)
	  {
	    my @aa = split /\s+/, $1;
	    if ($option eq 1) {
		for ($i = 0; $i < scalar @aa; $i ++) { $aa[$i] ++;} # adjust to get the edge number
	    }

	    my %bb = @aa;
	    @{$graph_obj{outedges}{$vertex_no}} = keys %bb;
	    for ($i = 0; $i < scalar @aa; $i += 2)
	      {
		# do isomorphism assignment of edges;
		$watson = $aa[$i];
		$crick  = $aa[$i+1];
#		$watson = $aa[$i] +1 if $option eq 2;
#		$crick  = $aa[$i+1]+1 if $option eq 2;
		$graph_obj{pal}{$watson} = $crick;
		$graph_obj{pal}{$crick} = $watson;
		
		$graph_obj{edge}{$watson}{from} = $vertex_no;
	      }
	  }
      }
    
    foreach my $e (keys %{$graph_obj{pal}})
    {
	$graph_obj{node_pal}{$graph_obj{edge}{$e}{from}} = $graph_obj{edge}{$graph_obj{pal}{$e}}{to};
	$graph_obj{node_pal}{$graph_obj{edge}{$e}{to}} = $graph_obj{edge}{$graph_obj{pal}{$e}}{from};
    }

    if ($verbose) {
      print "Graph Loaded:";
      printf "\t %d nodes; ", (scalar keys %{$graph_obj{node}});
      printf " %d edges. \n", (scalar keys %{$graph_obj{edge}});
    }
    return %graph_obj;
  }



#################################################
# sub WriteDot
# create .dot file from %graph, from rep_graph output
#################################################
sub WriteDot
{
    my ($graph_, $file) = @_;
    open (OUT, ">$file") or die;

    my @color = qw (black black red blue green);

    print OUT "digraph G {\n";
    print OUT "\tsize=\"8,8\";\n";

    foreach my $e (keys %{$$graph_{edge}})
    {
	my $style = ""; my $color_string = ""; my $color;
	$style = "style=bold," if ($$graph_{'len'}{$e} < 1000) ;
	if ( (my $i = $$graph_{'mul'}{$e}) > 1)
	{
	    if ( $i <= 3)
	    {
		$color = $color[$i];
	    }
	    else
	    {
		$color = 'green';
	    }
	    $color_string = "color=$color," ;
	}
	printf OUT "\t%d -> %d [$color_string $style label = \"%d,%d(%d)\"];\n",
	$$graph_{edge}{$e}{from}, $$graph_{edge}{$e}{to},
	$$graph_{'len'}{$e}, $e, $$graph_{'mul'}{$e};
    }
    print OUT "}\n";
    close OUT;

}

#################################################
# sub LoadDot 
# generate %graph, %len, %mul from rep_graph output
#################################################
sub LoadDot
{
    my ($dot_file, $verbose, $debug) = @_;
    my %found;
    my ($from, $to);
    my %graph;

    open (DOT, $dot_file) or die;
    while (<DOT>)
    {
	if (/(\d+) -> (\d+)/)
	{
	    $from = $1;
	    $to = $2;
	    $graph{node}{$from} = 1;
	    $graph{node}{$to} = 1;
	}

	if (/label = \"(\d+),(\d+)\((\d+)\)\"/)
	{
	    $graph{'mul'}{$2} = $3;
	    $graph{'len'}{$2} = $1;
	    $found{$2} = 1 if ($graph{'mul'}{$2} eq 1); 
	    $graph{edge}{$2}{from} = $from;
	    $graph{edge}{$2}{to} = $to;
	    push @{$graph{inedges}{$to}}, $2;
	    push @{$graph{outedges}{$from}}, $2;
	}
    }
    printf "graph loaded from dot file $dot_file, %d nodes, %d edges\n",
    scalar keys %{$graph{node}}, scalar keys %{$graph{edge}} if ($verbose);
    return %graph;
}



################################################
# sub LoadPath
################################################

sub LoadPath
{
    my $pfile = shift; # txt file with node and edge are listed, separated by space
    my %path;
    my $path_string= "";
    my ($e, $n);

    open (PATH, $pfile) or die;
    while (<PATH>)
    {
	chomp;
	$path_string .= " ".$_;
    }
    close PATH;
    
    my @crude = split /\s+/, $path_string;
    $n = shift @crude;
    push  @{$path{'node'}}, $n;
    my $e_index = 0;
    while ($e = shift @crude)
    {
	push @{$path{'edge'}}, $e;
	$path{'edge_index'}{$e} = $e_index++;
	$n = shift @crude;
	push @{$path{'node'}}, $n;
    }
    return %path;
}
################################################
# sub LoadPathGraph
################################################

sub LoadPathGraph
{
    my ($pfile, $pal_, $debug, $verbose, $option) = @_; 
    # $pfile = txt file with node and edge are listed, separated by space
    my %path;
    my $path_string= "";
    my ($e, $n) = (0, undef);

    open (PATH, $pfile) or die;
    while (<PATH>)
    {
	chomp;
	$path_string .= " ".$_;
    }
    close PATH;
    
    my @crude = split /\s+/, $path_string;
    while (undef eq $n) {$n = shift @crude};
    push  @{$path{'node'}}, $n;
    my $e_index = 0;
    while (($e = shift @crude) ne undef)
    {
      $e ++ if $option;
	push @{$path{'edge'}}, $e;
	$path{'edge_index'}{$e} = $e_index++;
	$n = shift @crude;
	push @{$path{'node'}}, $n;
    }
    
	push @{$path{'edge'}}, -1;
	$path{'edge_index'}{-1} = $e_index++;
    
    my $len = scalar @{$path{'edge'}} - 1;
    for (my $i = $len-1; $i >= 0; $i --)
    {
	my $e = $path{'edge'}[$i];
	push @{$path{'edge'}}, $$pal_{$e};
	$path{'edge_index'}{$$pal_{$e}} = $e_index++;
    }

    printf "total path loaded from $pfile: including %d edges\n", scalar @{$path{'edge'}} if $verbose;
    return %path;
}

################################################
# sub LoadReadPath
################################################

sub LoadReadPath
{
    my %path;
    my ($path_file) = @_;
    my ($ri, $l1, $l2, $num_edges, $dir, $to_read_next_line);
    my @edges;
    my @es;
    my %eh;

    open (PATH, $path_file) or die;
    while (<PATH>)
    {
	if (/SEQPATH\s(\d+)\s(\d+)\s(\d+)\s(\d+)/)
	{
	    $ri = $1;
	    ($l1, $l2) = ($2, $3);
	    $num_edges = $4;
	    $dir = 'f';
	    $to_read_next_line = 1;
	    next;
	}

	if (/REVPATH\s(\d+)\s(\d+)\s(\d+)\s(\d+)/)
	{
	    $ri = $1;
	    ($l1, $l2) = ($2, $3);
	    $num_edges = $4;
	    $dir = 'r';
	    $to_read_next_line = 1;
	    next;
	}

	if ($to_read_next_line)
	{
	    %eh = split; 
	    @es =  keys %eh;
	    if (scalar @es ne $num_edges)
	    {
		print "wrong path file format in $path_file: ", join ",", @es, "\n";
	    }
	    else
	    {
		@{$path{$ri}{$dir}} =  @es;
	    }
	    $to_read_next_line = 0;
	}

    }
    return %path;
}
 


################################################
# sub SaveReads
################################################
sub SaveReads
{
	my ($read_hash_, $read_list_, $outfile) = @_;
	my $out = Bio::SeqIO->new('-file' => ">$outfile",
				 										'-format' => 'Fasta');
	foreach my $rid (keys %$read_list_)
	{
                next unless $$read_hash_{$rid};
		$out->write_seq($$read_hash_{$rid});
	}
}


################################################
# sub LoadReads
################################################

sub LoadReads
  {
    use FileHandle;
    my ( $reads_file,$verbose,$debug ) = @_;
      my $step = 1000;
    print "Loading reads from file $reads_file (.=$step reads)\n" if $verbose;

    my $reads_io  = Bio::SeqIO->new('-file' => "$reads_file",
				 '-format' => 'Fasta');
    my %read_hash;

    my $i = 0;
    while ( my $readseq = $reads_io->next_seq() ) {
      $read_hash{'rid'}{$readseq->id} = $readseq;
      $read_hash{'array'}[$i] = $readseq->id;
      $read_hash{'index'}{$readseq->id} = $i;
      $read_hash{'length'}{$readseq->id} = $readseq->length; 
      $i ++;
      if (($verbose) and not ($i % $step))
      {
	  print ".";
	  STDOUT->autoflush(1);
      }
    }
    print "\n# reads loaded:", scalar (keys %{$read_hash{'rid'}}), "\n\n" if $verbose;
   
 return %read_hash;
  }

################################################
# sub LoadReadsWithQuals
################################################

sub LoadReadsWithQuals
  {
    use FileHandle;
    my ( $reads_file,$verbose,$debug,$quality_file,$id_field ) = @_;
      my $step = 1000;
    print "Loading reads from file $reads_file (.=$step reads)\n" if $verbose;

    my $reads_io  = Bio::SeqIO->new('-file' => "$reads_file",
				 '-format' => 'Fasta');
    my $quals_io  =  Bio::SeqIO->new('-file' => "$quality_file",
				 '-format' => 'qual');
    my %read_hash;

    my $i = 0;
    while ( my $readseq = $reads_io->next_seq() ) {
      my $readqual = $quals_io->next_seq();
  #    assert($readseq->id eq $readqual->id);
      my $id;
      if ($id_field eq 'id')
	{
	  $id=$readseq->id;
	}
      elsif ($id_field eq 'desc')
	{
	  $id=$readseq->desc;
	}
      else
	{
	  die "wrong id field: $id_field in loading reads\n";
	}
      $read_hash{'rid'}{$id} = $readseq;
      $read_hash{'qual'}{$id} = $readqual;
      $read_hash{'array'}[$i] = $id;
      $read_hash{'index'}{$id} = $i;
      $read_hash{'length'}{$id} = $readseq->length; 
      $i ++;
      if (($verbose) and not ($i % $step))
      {
	  print ".";
	  STDOUT->autoflush(1);
      }
    }
    print "\n# reads loaded:", scalar (keys %{$read_hash{'rid'}}), "\n\n" if $verbose;
   
 return %read_hash;
  }

################################################
# sub LoadGenome
################################################

sub LoadGenome
  {
    my ($genome_file, $option, $debug) = @_;
    my $genome_io  = Bio::SeqIO->new('-file' => "$genome_file",
				'-format' => 'Fasta');
    my $genome_seq = $genome_io->next_seq();

    if ($option eq 'm')
      {
	while (my $seq =  $genome_io->next_seq())
	  {
	    $genome_seq->seq($genome_seq->seq . $seq->seq);
	  }
      }

    printf "genome loaded: length = %d\n", $genome_seq->length if $debug;
return $genome_seq;
  }

################################################
# sub LoadSequence
################################################

sub LoadSingleSequence
  {
    my $sequence_file = shift;
    my $seq_io  = Bio::SeqIO->new('-file' => "$sequence_file",
				'-format' => 'Fasta');
    $seq_io->next_seq();
  }


################################################
# sub LoadThread
################################################

sub LoadThread
  {
    my ($thread_file, %edge_hash) =@_;
    my $forward_edge;
    my $reverse_edge;
my $edge_name;
my $euler_length;
my (%multiplicity, %start, %stop, %len, %strand);
my @tangle_edges;
    open (THREAD, $thread_file) or die "can not open thread file\n";
    while(<THREAD>)
      {
	if (/\(\s*(\d+)-\s*(\d+)\)\s*edge(\d+)\s*edge(\d+)\s*(\d+)/)
	  {
	    $forward_edge = $3;
	    $reverse_edge = $4;
	    $edge_name = 'edge'.$forward_edge;
	    $euler_length = $edge_hash{$edge_name}->length;

	    if ($5 == $euler_length) 	# authentize edge match by length
	      {
		$multiplicity{$forward_edge} ++;
		$start{$forward_edge}{$multiplicity{$forward_edge}} = $1;
		$stop{$forward_edge}{$multiplicity{$forward_edge}}  = $2;
		$len{$forward_edge}{$multiplicity{$forward_edge}} = $5;
		$strand{$forward_edge}{$multiplicity{$forward_edge}} = 'forward';
		
		$multiplicity{$reverse_edge} ++;
		$start{$reverse_edge}{$multiplicity{$reverse_edge}} = $1;
		$stop{$reverse_edge}{$multiplicity{$reverse_edge}}  = $2;
		$len{$reverse_edge}{$multiplicity{$reverse_edge}} = $5;
		$strand{$reverse_edge}{$multiplicity{$reverse_edge}} = 'reverse';
	      }
	  }
      }
    close THREAD;
    
    foreach my $edge (keys %multiplicity)
      {
	push @tangle_edges, $edge if $multiplicity{$edge} >= 2;
      }

  }


sub LoadReadCoordinate
{
    my ($read_coordinate_file, $read_ends_, $read_coordinate_ ) = @_;
    my $index;
    my $rid;
    my $read_pos;
    my $read_coord;
    my $read_count = 0;

    print "Loading read coordinates...\n";
    open (COOR, $read_coordinate_file) or die;
    while (<COOR>)
    {
	if (/^(.*): (\d+)\s*$/)
	{
	  $read_count ++;
	    $rid = $1;
	    $index = 0;
	}
	next unless $$read_ends_{$rid};

	if (/^(\d+):(\s+\d+\s-?1;.*)$/)
	{
	    $read_pos = $1;
	    my @coords = split ";", $2;
	    foreach $read_coord (@coords)
	    {
		if ($read_coord =~ /(\d+) (-?1)/)
		{
		    $$read_coordinate_{$rid}{'g_to_r'}{$1} = $read_pos;
		    $$read_coordinate_{$rid}{'r_to_g'}{$read_pos} = $1;
		    $$read_coordinate_{$rid}{'gstrand'}{$read_pos} = $2;
		    $index ++;
		}
	    }
	}
    }

    print scalar keys %$read_coordinate_, " reads has coordinates loaded\n";
}

sub LoadReadMapFile
{
    my ($read_map_file, $read_map_) = @_;
    my $rid;
    my $index;
    my $read_count = 0;

    open (RMAP, $read_map_file) or die;
    while (<RMAP>)
    {
	if (/^(.*): (\d+)/)
	{
	    $rid = $1;
	    $read_count ++;
	    $index = 0;
	}
	if (/^(-?1) g\((\d+),(\d+)\) r\((\d+),(\d+)\)/)
	{
	    $$read_map_{$rid}{'copy'}[$index]{'strand'} = $1;
	    $$read_map_{$rid}{'copy'}[$index]{'gstart'} = $1>0?$2:$3;
	    $$read_map_{$rid}{'copy'}[$index]{'gend'} = $1>0?$3:$2;
	    $$read_map_{$rid}{'copy'}[$index]{'rstart'} = $4;
	    $$read_map_{$rid}{'copy'}[$index]{'rend'} = $5;
	    $index ++;
	}
    }
    close RMAP;
    print "$read_count reads has been scanned\n";
    print scalar keys %$read_map_, " reads has been mapped \n";
}



sub LoadEdgeInterval
{
    my ($file, $intv_, $edge_map_, $total_len, $verbose) = @_;
    open (EI, $file) or die;
    my $e;
    my $ei;
  
    while (<EI>)
    {
	if (/EDGE (\d+) Length (\d+) Multiplicity (\d+)\./)
	{
	    $e = $1;
	    $ei = 0;
	    $$edge_map_{$e}{'length'} = $2;
	    $$edge_map_{$e}{'mul'} = $3;
	}

	if (/INTV (\d+) (\d+) (\d+) 0/)
	{
# 	  my $strand = 1-2*$1;
	  my $strand = $1;
	  $$edge_map_{$e}{'intv'}[$ei]{'strand'} = $strand;
	  my %thing;
	  if ($strand >= 0)
	    {
	      $$edge_map_{$e}{'intv'}[$ei]{'start'} = $2;
	      $$edge_map_{$e}{'intv'}[$ei]{'end'} = $2+$3;

	      $thing{'edge'} = $e;
	      $thing{'strand'} = $strand;
	      $thing{'start'} = $2;
	      $thing{'end'} = $2+$3;
	    }
	  else
	    {
	      $$edge_map_{$e}{'intv'}[$ei]{'start'} = $total_len-$2-1;
	      $$edge_map_{$e}{'intv'}[$ei]{'end'} = $total_len-$2-1-$3;

	      $thing{'edge'} = $e;
	      $thing{'strand'} = $strand;
	      $thing{'start'} = $total_len-$2-1-$3;
	      $thing{'end'} = $total_len-$2-1;
	    }
	  push @$intv_, \%thing;
	  $ei ++;
	}
    }

    if ($verbose)
      {
	print scalar keys %$edge_map_, " edges has been mapped \n";
	print scalar @$intv_, " intervals has been mapped \n";
      }
}

sub LoadEdgeInterval_Hash
{
    my ($file, $intv_, $edge_map_, $hash_, $total_len, $verbose) = @_;
    open (EI, $file) or die;
    my $e;
    my $ei;
  
    while (<EI>)
    {
	if (/EDGE (\d+) Length (\d+) Multiplicity (\d+)\./)
	{
	    $e = $1;
	    $ei = 0;
	    $$edge_map_{$e}{'length'} = $2;
	    $$edge_map_{$e}{'mul'} = $3;
	}

	if (/INTV (\d+) (\d+) (\d+) 0/)
	{
# 	  my $strand = 1-2*$1;
	  my $strand = $1;
	  $$edge_map_{$e}{'intv'}[$ei]{'strand'} = $strand;
	  my %thing;
	  if ($strand >= 0)
	    {
	      $$edge_map_{$e}{'intv'}[$ei]{'start'} = $2;
	      $$edge_map_{$e}{'intv'}[$ei]{'end'} = $2+$3;

	      $thing{'edge'} = $e;
	      $thing{'strand'} = $strand;
	      $thing{'start'} = $2;
	      $thing{'end'} = $2+$3;
	    }
	  else
	    {
	      $$edge_map_{$e}{'intv'}[$ei]{'start'} = $total_len-$2-1;
	      $$edge_map_{$e}{'intv'}[$ei]{'end'} = $total_len-$2-1-$3;

	      $thing{'edge'} = $e;
	      $thing{'strand'} = $strand;
	      $thing{'start'} = $total_len-$2-1-$3;
	      $thing{'end'} = $total_len-$2-1;
	    }
	  push @$intv_, \%thing;

	  $$hash_{strand}{$strand}{edge}{$e} = 1;
	  $ei ++;
	}
    }

    if ($verbose)
      {
	print scalar keys %$edge_map_, " edges has been mapped \n";
	print scalar @$intv_, " intervals has been mapped \n";
      }
}



1;

__END__
