***************
COMREP Package
***************

This package contains perl scripts for preparing, generating, and analyzing repeat domain graphs.

Procedures
1. Input: repeat.fa, a repeat library
2. Run ABA (1.01+) over the input repeat library, 
   - using the alignment program of your choice (blastn or cross_match) and the corresponding parser to generate glues for ABA.
3. Isolate connected components
   - supergraph.pl
   - subgraph.pl
4. Move subg*.dot to a subdirectory, and collect graph statistics
   - count_components.pl
   - count_singleton.pl
5. Simplify the subgraphs
   - beautify.pl, which invoke beautify_graph.pl
6. Discover Y-fork structures (in the beautified graphs)
   - Y_fork_finder.pl
7. Discover domains shared between repeat families of different biological orgins.
   - bio_CR_finder.pl, the biological origin of a repeat family is derived from the first two letters of its name. 


Degui Zhi
2005.07.01