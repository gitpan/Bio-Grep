#!perl -T 
use strict;
use warnings;
use lib 't';
use BioGrepTest;

use Test::More; 

my %prereq = BioGrepTest::check_prereq();
if (!$prereq{bioperl}) {
    plan skip_all => 'Bioperl not found';
}
elsif (!$prereq{bioperl_run}) {
    plan skip_all => 'Bioperl-run not found';
}

my $backendname  = 'Vmatch';
plan skip_all => 'Vmatch binary not in path' if
    BioGrepTest::find_binary_in_path( lc($backendname) ) eq '';

# make taint happy    
BioGrepTest::set_path( ( 'vmatch' ) );

plan tests => 2;

my $code =<<'EOT'
  use Bio::Grep;
  
  my $search_obj = Bio::Grep->new('Vmatch');	
  
  # $sbe is now a reference to the back-end
  my $sbe = $search_obj->backend;
 
  # define the location of the suffix arrays
  $sbe->settings->datapath('data');
  
  mkdir($sbe->settings->datapath);	
  
  # now generate a suffix array. you have to do this only once.
  $sbe->generate_database('t/Test.fasta', 'Description for the test Fastafile');
  
  # search in this suffix array
  $sbe->settings->database('Test.fasta');
  
  # search for the reverse complement and allow 2 mismatches
  $sbe->settings->query('UGAACAGAAAG');
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(2);

  # or you can use Fasta file with queries
  # $sbe->settings->query_file('Oligos.fasta');

  # $sbe->search();

  # Alternatively, you can specify the settings in the search call.
  # This also resets everything except the paths and the database
  # (because it is likely that they don't change when search is called
  # multiple times)

  $sbe->search( { query  =>  'UGAACAGAAAG',
                  reverse_complement => 1,
                  mismatches         => 2,
                 });  
  
  my @ids;

  # output some informations! 
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->alignment_string() . "\n\n";
     push @ids, $res->sequence_id;
  }
  
  # get the gene sequences of all matches as Bio::SeqIO object.
  # (to generate a Fasta file for example)
  my $seqio = $sbe->get_sequences(\@ids);

EOT
;

eval $code;
ok(!$@,"SYNOPSIS compiles") || diag $@;

$code = 'bllll';
eval $code;
ok($@,"bllll not compiles");


1;

# vim: ft=perl sw=4 ts=4 expandtab
