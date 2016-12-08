#this is an example of making a fasta database for quick indexing 

use Bio::Index::Fasta; # using fasta file format
use strict; # some users have reported that this is necessary

 
my $Index_File_Name = "$ARGV[0]";
 
my $inx = Bio::Index::Fasta->new(
    -filename => $Index_File_Name,
    -write_flag => 1);
 
$inx->make_index("$ARGV[1]")
