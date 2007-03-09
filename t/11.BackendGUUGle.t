#!perl -T 
################################################################################
# some backend tests, extra because GUUGle does not support mismatches
#
# Test fasta are sequences from ATH1.cdna, with MODIFIED sequences
################################################################################

use strict;
use warnings;

BEGIN{
    use lib 't';
    use Test::More; 
    use BioGrepTest;


    # make taint happy    
    BioGrepTest::set_path( ( 'guugle' ) );

    my %prereq = BioGrepTest::check_prereq();
    if (!$prereq{bioperl}) {
        plan skip_all => 'Bioperl not found';
    }
    elsif (!$prereq{bioperl_run}) {
        plan skip_all => 'Bioperl-run not found';
    }
}
my $backendname  = 'GUUGle';
plan skip_all => 'GUUGle binary not in path' if
BioGrepTest::find_binary_in_path( lc($backendname) ) eq '';
plan tests => 35;

use English qw( -no_match_vars );
use Cwd;
use Data::Dumper;
use Bio::Grep;

use Bio::Perl;

my %test_seq = (
    id   => 'At2g42200',
    desc => '68409.m05466 squamosa-promoter binding protein -related',
    seq  =>
        'accactctcgtctctttcttttttccttctgttctgtttctctctctaaacccaaaacagtcaaaatcagggaagccgaaattttctttgctttcttctcctttggtcctttctttaaacccgagacagttaggtttgtgtgagagagagaatgatgagtaaaaccctttctgtctgagtaagaggaaaccaacATGGAGATGGGTTCCAACTCGGGTCCGGGTCATGGTCCGGGTCAGGCAGAGTCGGGTGGTTCCTCCACTGAGTCATCCTCTTTCAGTGGAGGGCTCATGTTTGGCCAGAAGATCTACTTCGAGGACGGTGGTGGTGGATCCGGGTCTTCTTCCTCAGGTGGTCGTTCAAACAGACGTGTCCGTGGAGGCGGGTCGGGTCAGTCGGGTCAGATACCAAGGTGCCAAGTGGAAGGTTGTGGGATGGATCTAACCAATGCAAAAGGTTATTACTCGAGACACCGAGTTTGTGGAGTGCACTCTAAAACACCTAAAGTCACTGTGGCTGGTATCGAACAGAGGTTTTGTCAACAGTGCAGCAGGTTTCATCAGCTTCCGGAATTTGACCTAGAGAAAAGGAGTTGCCGCAGGAGACTCGCTGGTCATAATGAGCGACGAAGGAAGCCACAGCCTGCGTCTCTCTCTGTGTTAGCTTCTCGTTACGGGAGGATCGCACCTTCGCTTTACGAAAATGGTGATGCTGGAATGAATGGAAGCTTTCTTGGGAACCAAGAGATAGGATGGCCAAGTTCAAGAACATTGGATACAAGAGTGATGAGGCGGCCAGTGTCGTCACCGTCATGGCAGATCAATCCAATGAATGTATTTAGTCAAGGTTCAGTTGGTGGAGGAGGGACAAGCTTCTCATCTCCAGAGATTATGGACACTAAACTAGAGAGCTACAAGGGAATTGGCGACTCAAACTGTGCTCTCTCTCTTCTGTCAAATC'
);

my %hits_sequences = (
 'At1g01010.1' => 'ugaaaauggaggaucaaguuggguuug',
 'At1g01010.2' => 'aaaauggaggaucaaguuggguuug',
 'At1g01010.3' => 'ugaaaauggaggaucaaguuggguu',
);
my %hits_sequences2 = (
 'At1g01010.1' => 'aaauggaggaucaaguuggguuug',
 'At1g01010.2' => 'aaauggaggaucaaguuggguuug',
 'At1g01010.3' => 'aaauggaggaucaaguuggguu',
);

my $search_obj = Bio::Grep->new($backendname);
my $sbe        = $search_obj->backend;

# define own tmppath, so that we can check if all temporary files are deleted
$sbe->settings->tmppath('t/tmp');
mkdir("t/tmp");
mkdir("t/data");
BioGrepTest::delete_files;


$sbe->settings->reverse_complement(0);
$sbe->settings->datapath('t/data');
$sbe->generate_database_out_of_fastafile( 't/Test.fasta',
 'Description for Test.fasta' );
$sbe->generate_database_out_of_fastafile( 't/TestGUUGleExtend.fasta',
 'Description for Test.fasta' );
$sbe->settings->query('auggaggaucaaguugg');
$sbe->settings->database('TestGUUGleExtend.fasta');
$sbe->settings->gumismatches(0);
$sbe->search();
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, $sbe->settings->query, 'subject is query'); 
}
# upstream downstream tests
$sbe->settings->upstream(5);
$sbe->settings->downstream(5);
$sbe->search();
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, $sbe->settings->query, 'subject is query'); 
    is($res->sequence->seq, $hits_sequences{$res->sequence->id}, 'sequence correct'); 
}
$sbe->settings->query_length(14);
$sbe->search();
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, $sbe->settings->query, 'subject is query'); 
    is($res->sequence->seq, $hits_sequences{$res->sequence->id}, 'sequence correct'); 
}
$sbe->settings->query_length(14);
$sbe->settings->query('cccgaggaucaaguugg');
$sbe->search();
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, 'gaggaucaaguugg', 'subject is query'); 
    is($res->sequence->seq, $hits_sequences2{$res->sequence->id}, 'sequence correct'); 
}
$sbe->settings->query('cccgaggaucaaguuuu');
$sbe->search();
while (my $res = $sbe->next_res() ) {
                           #ccaacuugauccuc  # rev_com
                           #cuccuaguucaacc  # complement
    is($res->subject->seq, 'gaggaucaaguugg', 'subject is query'); 
    is($res->sequence->seq, $hits_sequences2{$res->sequence->id}, 'sequence correct'); 
}
# test reverse complement
$sbe->settings->reverse_complement(1);
my $query = 'auggaggaucaaguugg';
$query =~ s/u/t/g;
$sbe->settings->query(revcom_as_string($query));
$sbe->search();
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, 'auggaggaucaaguugg', 'subject is query'); 
}



# test database
$sbe->search({
    reverse_complement => 0,
    database           => 'Test.fasta',
    gumismatches => 0,
    query              => 'CAGAGTCGGGTGGTTCCTCCACTGAGTCATCCTCTTTCAGTGGAGGGCTCAT',
});

my $test_seq_internal_id;
while (my $res = $sbe->next_res() ) {
     $test_seq_internal_id = $res->sequence_id
}
isnt( $test_seq_internal_id, '', 'Found internal id' ) ;
my $seqio        = $sbe->get_sequences( [$test_seq_internal_id] );
my $test_seq_obj = $seqio->next_seq();

SKIP: {
    skip "Could not get sequence object ($backendname)", 3
        if !defined $test_seq_obj;

    is( $test_seq_obj->id,   $test_seq{id} );
    is( $test_seq_obj->desc, $test_seq{desc} );
    is( $test_seq_obj->seq,  $test_seq{seq} );
}

# test for GUUGle specific exceptions
$sbe->settings->upstream(10);
$sbe->settings->downstream(5);
$sbe->settings->gumismatches(0);

eval { $sbe->search(); };

ok($EVAL_ERROR, 'Exception occured with different values for up- and ' .
     'downstream.');
### 

$sbe->settings->downstream(10);

eval { $sbe->search(); };

ok(!$EVAL_ERROR, 'No exception occured with equal values for up- and ' .
     'downstream.') || diag $EVAL_ERROR;
###

$sbe->settings->upstream_reset;

eval { $sbe->search(); };

ok($EVAL_ERROR, 'Exception occured with undef up- and ' .
     'def. downstream.');

###
 
eval { $sbe->search( { query_file => 't/Test.fasta',  
                 query_length => 20,
                 gumismatches => 0,
             } ); 
     }; 

ok($EVAL_ERROR, 'Exception occured when revcom not set');

###

eval { $sbe->search( { query_file => 't/Test.fasta', 
                 query_length => 20, 
                 gumismatches => 0,
                 reverse_complement => 1  } 
             ); 
     }; 

ok(!$EVAL_ERROR, 'No exception occured when revcom set') || diag $EVAL_ERROR;

###

diag "\nNow you should see a warning\n";

eval { $sbe->search( { query_file => 't/Test.fasta', 
                 query_length => 20, 
                 gumismatches => 1,
                 reverse_complement => 1  } 
             ); 
     }; 

ok(!$EVAL_ERROR, 'No exception occured when revcom set') || diag $EVAL_ERROR;

###

eval { $sbe->search( { query_file => 't/Test.fasta', 
                       gumismatches => 0,
                       reverse_complement => 1  } 
             ); 
     }; 

ok($EVAL_ERROR, 'Exception occured when query_length is missing');
             
1;

# vim: ft=perl sw=4 ts=4 expandtab
