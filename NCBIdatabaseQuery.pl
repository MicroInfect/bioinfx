#!/usr/local/bin/perl 
use Bio::DB::EUtilities;
use strict;
use warnings;

my $line = ();
my $filename = 'names.txt';
open(my $fh ,'<' ,$filename) or die "could not open file '$filename' $!";


foreach $line (<$fh>)
        {
        chomp $line;
        my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
                                       -db     => 'nuccore',
                                       -term   => "$line AND srcdb_refseq[prop]",
                                       -email  => 'mymail@foo.bar',
                                       -retmax => 5000);

# query terms are mapped; what's the actual query?
print "Query translation: ",$factory->get_query_translation,"\n";

# query hits
print "Count = ",$factory->get_count,"\n";

# UIDs
my @ids = $factory->get_ids;
print "@ids";
}
print "\n";
