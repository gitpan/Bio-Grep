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


    my %prereq = BioGrepTest::check_prereq();
    if (!$prereq{bioperl}) {
        plan skip_all => 'Bioperl not found';
    }
    elsif (!$prereq{bioperl_run}) {
        plan skip_all => 'Bioperl-run not found';
    }
}
my $backendname  = 'RE';
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
 'At1g01010.1' => 'tgaaaatggaggatcaagttgggtttg',
 'At1g01010.2' => 'aaaatggaggatcaagttgggtttg',
 'At1g01010.3' => 'tgaaaatggaggatcaagttgggtt',
);
my %hits_sequences3 = (
 'At1g01010.1' => 'tgaaaatggaggatcaagttgggt',
 'At1g01010.2' => 'aaaatggaggatcaagttgggt',
 'At1g01010.3' => 'tgaaaatggaggatcaagttgggt',
);
my %hits_sequences4 = (
 'At1g01010.1' => 'aaatggaggatcaagttgggtttg',
 'At1g01010.2' => 'aaatggaggatcaagttgggtttg',
 'At1g01010.3' => 'aaatggaggatcaagttgggtt',
);
my %hits_sequences5 = (
 'At1g01010.1' => 'atggaggatcaagttgggtttg',
 'At1g01010.2' => 'atggaggatcaagttgggtttg',
 'At1g01010.3' => 'atggaggatcaagttgggtt',
);

my $sbe = Bio::Grep->new($backendname);

# define own tmppath, so that we can check if all temporary files are deleted
$sbe->settings->tmppath('t/tmp');
mkdir("t/tmp");
mkdir("t/data");
BioGrepTest::delete_files;


$sbe->settings->reverse_complement(0);
$sbe->settings->datapath('t/data');
$sbe->generate_database( 't/Test.fasta',
 'Description for Test.fasta' );
$sbe->generate_database( 't/TestGUUGleExtend.fasta',
 'Description for Test.fasta' );
$sbe->settings->query('ATGGAGGATCAAGTTGG');
$sbe->settings->database('TestGUUGleExtend.fasta');
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
    is(lc($res->sequence->seq), lc($hits_sequences{$res->sequence->id}), 'sequence correct'); 
}
# test reverse complement
$sbe->settings->reverse_complement(1);
my $query = 'ATGGAGGATCAAGTTGG';
$sbe->settings->query(revcom_as_string($query));
$sbe->search();
while (my $res = $sbe->next_res() ) {
    is(lc($res->subject->seq), lc($query), 'subject is query'); 
}

# different up/downstream
$sbe->search({
        query => 'ATGGAGGATCAAGTTGG',
        upstream => 5,
        downstream => 2,
    });
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, $sbe->settings->query, 'subject is query'); 
    is(lc($res->sequence->seq), lc($hits_sequences3{$res->sequence->id}), 'sequence correct'); 
}
$sbe->search({
        query => 'ATGGAGGATCAAGTTGG',
        upstream => 2,
        downstream => 5,
    });
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, $sbe->settings->query, 'subject is query'); 
    is(lc($res->sequence->seq), lc($hits_sequences4{$res->sequence->id}), 'sequence correct'); 
}
$sbe->search({
        query => 'ATGGAGGATCAAGTTGG',
        downstream => 5,
    });
while (my $res = $sbe->next_res() ) {
    is($res->subject->seq, $sbe->settings->query, 'subject is query'); 
    is(lc($res->sequence->seq), lc($hits_sequences5{$res->sequence->id}), 'sequence correct'); 
}

# test database
$sbe->search({
    reverse_complement => 0,
    database           => 'Test.fasta',
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

# test for RE specific exceptions

eval { $sbe->search( { query => 'GAA[AT]G', reverse_complement => 1 } ) };

cmp_ok($EVAL_ERROR, '=~', qr{Query does not look like a DNA/RNA sequence.}, 
    'Exception occured with revcom and regex') || diag $EVAL_ERROR;


1;

# vim: ft=perl sw=4 ts=4 expandtab
