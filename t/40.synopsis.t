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

plan tests => 5;

# VMATCH 

my $backendname  = 'Vmatch';

SKIP:{

    skip 'Vmatch binary not in path', 3 if
        BioGrepTest::find_binary_in_path( lc($backendname) ) eq '';

    # make taint happy    
    BioGrepTest::set_path( ( 'vmatch' ) );


my $code =<<'EOT'
  use Bio::Grep;
  
  my $search_obj = Bio::Grep->new('Vmatch');	
  
  # $sbe is now a reference to the back-end
  my $sbe = $search_obj->backend;
 
  # define the location of the suffix arrays
  $sbe->settings->datapath('t/data');
  
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

$code =<<'EOT'
  use Bio::Grep;
  use Bio::SeqIO;
 
  my $sbe = Bio::Grep->new('Vmatch')->backend;
 
  my $out = Bio::SeqIO->new( -format => 'Fasta',
                             -file   => '>motifs.fasta',
                           );
 
  # you have an array with DNA sequences
  my @motifs = ( 'aaaaaa', 'gggggg' );
  
  for my $i (0 .. $#motifs ) {
     my $seq = Bio::Seq->new(
             -id => $i,
             -seq => $motifs[$i],
         );
     $out->write_seq($seq);
  }
 
  $sbe->search({
     datapath   => 't/data',
     database   => 'Test.fasta',
     query_file => 'motifs.fasta',
     complete   => 1,
  });

EOT
;
    eval $code;
    ok(!$@,"Cookbook recipe motifs solution a compiles") || diag $@;

    unlink 'motifs.fasta';
}

# RE

my $code =<<'EOT'
  use Bio::Grep::Backend::RE;
  
  # configure our search back-end, in this case RE
  my $sbe = Bio::Grep::Backend::RE->new();
  
  $sbe->settings->datapath('t/data');
  
  # generate a database. you have to do this only once. 
  $sbe->generate_database('t/Test.fasta');
  
  $sbe->settings->database('Test.fasta');
  
  # search for the reverse complement 
  $sbe->settings->query('TGAACAGAAAGCTCATGAGCC');
  $sbe->settings->reverse_complement(1);
  
  $sbe->search();
  
  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";
  }

  # sequences with at least 10 As
  $sbe->search({ query => '[A]{10,}' });
 
  # some SNPs
  $sbe->search({query => '[CG]TGC[AT]CTCTTCT[CG]TCA'});


EOT
;
eval $code;
ok(!$@,"RE SYNOPSIS compiles") || diag $@;

$code =<<'EOT'
  use Bio::Grep;
 
  my $sbe = Bio::Grep->new('RE')->backend;
 
  my $motif = '[AC]{4}TAAAA[AGCT]GG';
 
  $sbe->search({
     datapath  => 't/data',
     database  => 'Test.fasta',
     query      => $motif,
  });
EOT
;
eval $code;
ok(!$@,"Cookbook recipe motifs solution b compiles") || diag $@;

$code = 'bllll';
eval $code;
ok($@,"bllll not compiles");

BioGrepTest::delete_files;

1;

# vim: ft=perl sw=4 ts=4 expandtab
