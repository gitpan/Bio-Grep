#!perl -T 
################################################################################
# some backend tests
#
# Test fasta are sequences from ATH1.cdna, with MODIFIED sequences
################################################################################

use strict;
use warnings;

BEGIN{
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
}

our %tests = ( Agrep => 64, Vmatch => 123, Hypa => 79, GUUGle => 17 );
my $number_tests = 1;

foreach (keys %tests) {
    $number_tests += $tests{$_};
}
plan tests => $number_tests;


# backends 
################################################################################


use English qw( -no_match_vars );
use Cwd;
use Scalar::Util qw/tainted/;
use Data::Dumper;



use Bio::Grep;
use Bio::Perl;

# the number of files we assume in the data directory after
# database generation.
my %backend_filecnt = ( Agrep => 6, Vmatch => 16, Hypa => 31, GUUGle => 4 );

# the hits in Test.fasta we assume with n mismatches
my %hits = (
    At1g53160 => 2,
    At2g42200 => 1,
    At1g27360 => 1,
    At3g47170 => 3,
    At1g22000 => 4,
    At1g09950 => 5
);

my %test_seq = (
    id   => 'At2g42200',
    desc => '68409.m05466 squamosa-promoter binding protein -related',
    seq  =>
        'accactctcgtctctttcttttttccttctgttctgtttctctctctaaacccaaaacagtcaaaatcagggaagccgaaattttctttgctttcttctcctttggtcctttctttaaacccgagacagttaggtttgtgtgagagagagaatgatgagtaaaaccctttctgtctgagtaagaggaaaccaacATGGAGATGGGTTCCAACTCGGGTCCGGGTCATGGTCCGGGTCAGGCAGAGTCGGGTGGTTCCTCCACTGAGTCATCCTCTTTCAGTGGAGGGCTCATGTTTGGCCAGAAGATCTACTTCGAGGACGGTGGTGGTGGATCCGGGTCTTCTTCCTCAGGTGGTCGTTCAAACAGACGTGTCCGTGGAGGCGGGTCGGGTCAGTCGGGTCAGATACCAAGGTGCCAAGTGGAAGGTTGTGGGATGGATCTAACCAATGCAAAAGGTTATTACTCGAGACACCGAGTTTGTGGAGTGCACTCTAAAACACCTAAAGTCACTGTGGCTGGTATCGAACAGAGGTTTTGTCAACAGTGCAGCAGGTTTCATCAGCTTCCGGAATTTGACCTAGAGAAAAGGAGTTGCCGCAGGAGACTCGCTGGTCATAATGAGCGACGAAGGAAGCCACAGCCTGCGTCTCTCTCTGTGTTAGCTTCTCGTTACGGGAGGATCGCACCTTCGCTTTACGAAAATGGTGATGCTGGAATGAATGGAAGCTTTCTTGGGAACCAAGAGATAGGATGGCCAAGTTCAAGAACATTGGATACAAGAGTGATGAGGCGGCCAGTGTCGTCACCGTCATGGCAGATCAATCCAATGAATGTATTTAGTCAAGGTTCAGTTGGTGGAGGAGGGACAAGCTTCTCATCTCCAGAGATTATGGACACTAAACTAGAGAGCTACAAGGGAATTGGCGACTCAAACTGTGCTCTCTCTCTTCTGTCAAATC'
);

my %hits_sequences = (
    At1g27360 => 'tgaatctcaagatatccaccGTGCTCTCTCTCTTCTGTCAacctcttcgg',
    At2g42200 => 'gggaattggcgactcaaactGTGCTCTCTCTCTTCTGTCAaatc'
    ,    # test for too short downstream
    At1g53160 => 'agattagatagagaagctgtCTGCTCTCTCTCTTCTGTCAtctaaacttc',
    At3g47170 => 'ttgataaaggcaacggcctcGTGCACTCTCTCTTCTCTCAattacttgga',
    At1g22000 => 'taccccgatgaaaagtttctGAGATCTCTTTCTTCTGTCAaacatctctt',
    At1g09950 =>
        'gccggggacaacGTTTTCACTTTCTTCTGCCCaccgtggttt'    # too short upstream
);
my $tests = 0;

my %sort_modes = (
    GUUGle  => [ 'ga', 'gd' ],
    Hypa    => [ 'ga', 'gd' ],
    Agrep   => [ ],
    Vmatch  => [ sort(              
               qw(la ld ia id ja jd ea ed sa sd ida idd ga gd) ) ],
);


BACKEND:
foreach my $backendname ( sort keys %tests ) {
SKIP: {
        if ( BioGrepTest::find_binary_in_path( lc($backendname) ) eq '' ) {
            skip "$backendname not found in path", $tests{$backendname};
        }
        else {
            $tests++;
        }
        diag("\n*** Testing $backendname ***");
        BioGrepTest::set_path( ( map { lc($_) } keys %backend_filecnt),
            'RNAcofold' );
        my $search_obj = Bio::Grep->new($backendname);
        my $sbe        = $search_obj->backend;
        my %asm        = $sbe->available_sort_modes();
        is_deeply([ sort keys %asm ], 
                  $sort_modes{$backendname}, 'sortmodes as expected');

 # define own tmppath, so that we can check if all temporary files are deleted
        $sbe->settings->tmppath('t/tmp');
        mkdir("t/tmp");
        mkdir("t/data");
        BioGrepTest::delete_files;
        

        $sbe->settings->datapath('t/data');
        eval {
        $sbe->generate_database_out_of_fastafile( 't/wrong\ named.fasta',
             'Description for wrong\ named.fasta' );
        };
        ok( $EVAL_ERROR, "exception occured with invalid filename" );
        
        eval {
            $sbe->generate_database_out_of_fastafile( 't/wrong.fasta',
                'Description for wrong.fasta' );
        };
        ok( $EVAL_ERROR, "exception occured with not existing file" );

        $sbe->settings->datapath('t/wrongdata');
        eval {
            $sbe->generate_database_out_of_fastafile( 't/Test.fasta',
                'Description for Test.fasta' );
        };
        ok( $EVAL_ERROR, "exception occured with not existing datapath" );

        $sbe->settings->datapath('t/data');

        eval {
            $sbe->generate_database_out_of_fastafile( 't/Test.fasta',
                'Description for Test.fasta' );
        };
        ok( !$EVAL_ERROR, "no exception occured with dna fasta ($backendname)" );

        is($sbe->get_alphabet_of_database('Test.fasta'), 'dna', "Test.fasta dna
        ($backendname)");

        # test if all files are there after database generation
################################################################################
        opendir( DIR, "t/data" ) || die "can't opendir t/data: $!";
        my @files = grep { /Test/ && -f "t/data/$_" } readdir(DIR);
        closedir DIR;

        foreach my $file (@files) {
            if ( $file =~ /nfo$/ ) {
                open( FILE, "t/data/$file" );
                my $filecontent = '';
                while (<FILE>) {
                    $filecontent .= $_;
                }
                close(FILE);
                is( $filecontent,
                    'Description for Test.fasta',
                    'Description ok'
                );
            }
        }
        is( scalar(@files), $backend_filecnt{$backendname},
            "$backend_filecnt{$backendname} $backendname index files created"
        );

        eval {
            $sbe->generate_database_out_of_fastafile( 't/Test_pep.fasta',
                'Description for Test_pep.fasta' );
        };
        
        if ($backendname eq 'GUUGle') {
            ok( $EVAL_ERROR, "exception occured with peptide fasta
                ($backendname)" );
        }
        else {
            ok( !$EVAL_ERROR, "no exception occured with peptide fasta
                ($backendname)" );
        }
        
        opendir( DIR, "t/data" ) || die "can't opendir t/data: $!";
        my @files2 = grep { /Test_pep/ && -f "t/data/$_" } readdir(DIR);
        closedir DIR;

        if ($backendname ne 'GUUGle') {
            foreach my $file (@files2) {
                if ( $file =~ /nfo$/ ) {
                    open( FILE, "t/data/$file" );
                    my $filecontent = '';
                    while (<FILE>) {
                        $filecontent .= $_;
                    }
                    close(FILE);
                    is( $filecontent,
                        'Description for Test_pep.fasta',
                        'Description ok'
                    );
                }
            }
            is( scalar(@files2), $backend_filecnt{$backendname},
                "$backend_filecnt{$backendname} $backendname index files created"
            );
        }
        push @files, @files2;
        
        if (defined $sbe->features->{PROTEINS}) {

            is($sbe->get_alphabet_of_database('Test_pep.fasta'), 'protein', "Test.fasta
            protein
            ($backendname)");
            
            is_deeply( { $sbe->get_databases },
                { 'Test.fasta' => 'Description for Test.fasta', 'Test_pep.fasta'
                => 'Description for Test_pep.fasta' } );
        }

       
        $sbe->settings->database('Test.fasta');

        my $sequence ='ttattagatataccaaaccagagaaaacaaatacat';
        
        if (defined $sbe->features->{NATIVE_ALIGNMENTS}) {
            my $msu =
            'aaaTTATTAGATATACCAAACCAGAGAAAACAAATACATaatcggagaaatacagattacagagagcga';

            if (defined $sbe->features->{GUMISMATCHES}) {
                $sbe->settings->gumismatches(0);
            }    

            $sbe->settings->query($sequence);
            $sbe->settings->upstream(30);
            $sbe->settings->downstream(30);
            $sbe->search();
            while (my $res = $sbe->next_res) {
                my $t_msu = $res->mark_subject_uppercase();
                $t_msu =~ tr/Uu/Tt/;
                is($t_msu, $msu, 'mark subject uppercase works');
                my $t_el = 3 + length($sequence)+30;
                is(length($t_msu),$t_el,
                    'length of sequence correct');
                is($res->begin, 3,'begin is 3');
                my $t_ss = $res->sequence->seq;
                $t_ss =~ tr/Uu/tt/;
                is(lc($t_ss), lc($msu), 'subject ok');
                my $seqio = $sbe->get_sequences([$res->sequence_id ]);
                my $db_seq = $seqio->next_seq;
                is(lc($db_seq->subseq(4,3+length($sequence))),
                   lc($sequence),'subject found in db');
                is(lc($db_seq->subseq(1,$t_el)),
                   lc($t_ss),'sequence found in db');
            }   
        $sbe->settings->set({});    
        }     
        if (defined $sbe->features->{MAXHITS}) {
            $sequence ='ttatt';

            $sbe->settings->query($sequence);
            if (defined $sbe->features->{GUMISMATCHES}) {
                $sbe->settings->gumismatches(0);
            }    
            $sbe->settings->query($sequence);
            $sbe->settings->maxhits(5);
            $sbe->search();
            my $cnt = 0;
            while (my $res = $sbe->next_res) {
                $cnt++;
            }    
            is ($cnt, 5, 'maxhits(5) returned 5 hits');

            $sbe->settings->set({});   
        }

        goto CLEANUP if $backendname eq 'GUUGle';

        my $test_seq_internal_id = '';
        $sequence             = 'tgacagaagagagtgagcac';

        # now search for 1 to 5 mismatches, test if reverse complement works
################################################################################
        for my $j ( 0 .. 1 ) {
            if ( $j == 0 ) {
                $sbe->settings->query($sequence);
                $sbe->settings->reverse_complement(1);
            }
            else {
                $sbe->settings->query( revcom_as_string($sequence) );
                $sbe->settings->reverse_complement(0);
            }
            for my $i ( 1 .. 5 ) {
                # test some vmatch flags
                # they should not change anything in the results
                if ($j && $i == 5 && defined( $sbe->features->{ONLINE} )) {
                    $sbe->settings->online(1);
                }    
                if (!$j && $i == 5 && defined( $sbe->features->{ONLINE} )) {
                    $sbe->settings->online(0);
                }    
                if ($j && $i == 4 && defined( $sbe->features->{QSPEEDUP} )) {
                    $sbe->settings->qspeedup(2);
                }    
                if ($j == 0 && $i == 4 && defined( $sbe->features->{QSPEEDUP} )) {
                    $sbe->settings->qspeedup(2);
                }    
                $sbe->settings->mismatches($i);

                $sbe->search();
                $sbe->settings->online_reset;
                $sbe->settings->qspeedup_reset;
                my @ids = ();
                while (my $res = $sbe->next_res ) {

                    #   print STDERR "IS: " . $res->sequence->id . "\n";
                    push ( @ids, $res->sequence->id );
                    $test_seq_internal_id = $res->sequence_id
                        if $res->sequence->id eq $test_seq{id};
                    if ($i == 1 && $j == 1) {
                        is(uc($res->query->seq), uc($sbe->settings->query), 
                            'query sequence in results');
                    }    
                }
                my @shouldbe = _get_ids_where_mm_smaller($i);
                foreach (@shouldbe) {

                    #    print "SHOULDBE: $_\n";
                }
                @ids = sort @ids;
                is_deeply( \@shouldbe, \@ids , "Results ok" );
                #warn Dumper(\@shouldbe, \@ids);
            }
        }

        # test if backend returns correct sequence

        isnt( $test_seq_internal_id, '', 'Found internal id' );
        my $seqio        = $sbe->get_sequences( [$test_seq_internal_id] );
        my $test_seq_obj = $seqio->next_seq();

    SKIP: {
            skip "Could not get sequence object ($backendname)", 3
                if !defined $test_seq_obj;

            is( $test_seq_obj->id,   $test_seq{id} );
            is( $test_seq_obj->desc, $test_seq{desc} );
            is( $test_seq_obj->seq,  $test_seq{seq} );
        }

        # testing, if backends find matches at the end of long sequences
################################################################################
        $sbe->settings->query('CGAGCTGATGCAAAGCTCGCGGGACTGA');
        $sbe->settings->reverse_complement(0);
        $sbe->settings->mismatches(0);
        $sbe->search();
        
        my @ids = ();
        while (my $res = $sbe->next_res  ) {
            push ( @ids, $res->sequence->id );
        }
        my @shouldbe = qw( At1g67120 );
        ok( compare_arrays( \@shouldbe, \@ids ), "Results ok" );
        
        if (defined $sbe->features->{COMPLETE}) {
            $sbe->settings->complete(1);
        }    
        $sbe->search();
        
        @ids = ();
         while (my $res = $sbe->next_res  ) {
            push ( @ids, $res->sequence->id );
        }
        ok( compare_arrays( \@shouldbe, \@ids ), "Results ok" );

        if ($backendname eq 'Vmatch') {
            for my $i ( 0 .. 2) {
                $sbe->settings->query_file('t/Test_query.fasta');
                if ($i == 2) {
                    $sbe->settings->showdesc(20);
                }
                if ($i == 1) {
                    $sbe->settings->query_file('t/Test_query_revcom.fasta');
                    $sbe->settings->reverse_complement(1);
                }

                $sbe->search();
                
                @ids = ();
                my %multi_query_result = (
                    'At1g01020.1' => {id => 'b', desc => 'descb', seq =>
                        'CGAGTGTGAACGCATGATTATTTTCATCGATTTAA'}, 
                    'At1g01030.1' => {id => 'c', desc => 'descc',
                    seq => 'gttttcttccgattctagggttttcatatttc'},
                    'At1g01010.1' => {id => 'a', desc => 'desca', seq =>
                    'TGTAGTGAGGGCTTTCGTGGTAAGATT' }
                ); 
                while (my $res = $sbe->next_res ) {
                    
                    is($res->query->id,
                        $multi_query_result{$res->sequence->id}->{id});
                    is($res->query->desc,
                        $multi_query_result{$res->sequence->id}->{desc});
                    if ($sbe->settings->reverse_complement) {
                        is($res->query->revcom->seq,
                            $multi_query_result{$res->sequence->id}->{seq});
                    }
                    else {
                        is($res->query->seq,
                            $multi_query_result{$res->sequence->id}->{seq});
                    }    
                }
                $sbe->settings->showdesc_reset;
                $sbe->settings->reverse_complement_reset;
            }
        }

        $sbe->settings->complete_reset();
        $sbe->settings->query_file_reset();
################################################################################
        # test upstream/downstream
################################################################################
        if ( $backendname eq 'Agrep' ) {
            ok( !( defined( $sbe->features->{UPSTREAM} ) ),
                "$backendname does not support upstream"
            );
            ok( !( defined( $sbe->features->{DOWNSTREAM} ) ),
                "$backendname does not support downstream"
            );
        }
        else {
            ok( defined( $sbe->features->{UPSTREAM} ),
                "$backendname does support upstream"
            );
            ok( defined( $sbe->features->{DOWNSTREAM} ),
                "$backendname does support downstream"
            );

            my $sequence = 'tgacagaagagagtgagcac';
            if( defined( $sbe->features->{COMPLETE} )) {
                $sbe->settings->complete(0);
            }    
            $sbe->settings->query($sequence);
            $sbe->settings->reverse_complement(1);
            $sbe->settings->mismatches(5);
            $sbe->settings->no_alignments(1);
            $sbe->settings->upstream(20);
            $sbe->settings->downstream(10);

            $sbe->search();
            $sbe->settings->complete_reset;
            my @ids = ();
            while (my $res = $sbe->next_res  ) {

                #         print STDERR "IS: " . $res->sequence->id . "\n";
                #       print STDERR $res->mark_subject_uppercase() . "\n";
                my $subject = $hits_sequences{ $res->sequence->id };
                $subject =~ s/[^AGCT]//g;
                is( $res->subject->id, $res->sequence->id,
                    "Subject id same as sequence id" );
                is( uc( $res->subject->seq ), $subject, "Subject correct" );
                is( $res->mark_subject_uppercase(),
                    $hits_sequences{ $res->sequence->id },
                    "sequence und marking subject correct."
                );
                push ( @ids, $res->sequence->id );
            }
            my @shouldbe = _get_ids_where_mm_smaller(5);
            foreach (@shouldbe) {

                #    print "SHOULDBE: $_\n";
            }
            @ids = sort @ids;
        SKIP: {
                skip "Bug in downstream/upstream in hypa", 1
                    if $backendname eq "Hypa";
                ok( compare_arrays( \@shouldbe, \@ids ), "Results ok" );
            }

        }

        # check exceptions with insecure sortmode
        $sbe->settings->sort('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong sort mode" );
        $sbe->settings->sort('gd');
        {
            no warnings;
            eval { $sbe->search() };
            if ( $backendname eq "Agrep" ) {
                ok( $EVAL_ERROR,
                    "$backendname: No exception occured with valid sort mode (Agrep exception)"
                );
            }
            else {
                ok( !$EVAL_ERROR,
                    "$backendname: No exception occured with valid sort mode (Agrep exception)"
                );
            }
        }

        $sbe->settings->sort_reset;

        #check exceptions with insecure database
        $sbe->settings->database('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong database" );
        $sbe->settings->database('Test.fasta');
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct database" );
        $sbe->settings->database_reset;
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with no database" );

        # check exceptions with insecure mismatches
        $sbe->settings->database('Test.fasta');
        $sbe->settings->mismatches('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong mismatches" );
        $sbe->settings->mismatches(0);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct mismatches" );

        # check exceptions with insecure insertions
        $sbe->settings->insertions('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong insertions" );
        $sbe->settings->insertions(0);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct insertions" );

        # check exceptions with insecure query_length
        $sbe->settings->query_length('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong query_length" );
        $sbe->settings->query_length(10);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct query_length" );
        $sbe->settings->query_length_reset;
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct query_length" );

        # check exceptions with insecure upstream
        $sbe->settings->upstream('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong upstream" );
        $sbe->settings->upstream(0);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct upstream" );
        $sbe->settings->upstream_reset;
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct upstream" );

        # check exceptions with insecure downstream
        $sbe->settings->downstream('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong downstream" );
        $sbe->settings->downstream(0);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct downstream" );
        $sbe->settings->downstream_reset;
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct upstream" );

        # check exceptions with insecure maxhits
        $sbe->settings->maxhits('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong maxhits" );

        $sbe->settings->mismatches(5);
        $sbe->settings->maxhits(3);
        $sbe->settings->query($sequence);
        $sbe->settings->reverse_complement(1);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct maxhits" );
        if ( defined $sbe->features->{MAXHITS} ) {
            is( scalar( @{$sbe->results}  ),
                3, "only 3 results ($backendname)" );
        }
        $sbe->settings->maxhits_reset;
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct maxhits" );
        is( scalar(  @{$sbe->results}  ), 6, "all results ($backendname)" );
        $sbe->settings->mismatches(0);

        # check exceptions with insecure edit distance
        $sbe->settings->editdistance('&& ls *;');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured with wrong editdistance" );
        $sbe->settings->editdistance(1);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct editdistance" );
        $sbe->settings->editdistance_reset;
        $sbe->settings->mismatches(1);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured with correct editdistance" );
        {
            $sbe->settings->editdistance(1);
            $sbe->settings->mismatches(1);
            $sbe->verbose(2);
            eval { $sbe->search() };
            ok( $EVAL_ERROR, "warning occured ($backendname)" );
            $sbe->verbose(1);
        }
        $sbe->settings->database('Test_pep.fasta');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured when searching dna seq in pep db" );
        
        $sbe->settings->query('IDSPALKELHLSEV');
        $sbe->settings->reverse_complement(0);
        eval { $sbe->search() };
        ok( !$EVAL_ERROR, "No Exception occured when searching pep seq in pep db" );
        is(scalar @{$sbe->results},1);

        $sbe->settings->database('Test.fasta');
        eval { $sbe->search() };
        ok( $EVAL_ERROR, "Exception occured when searching pep seq in dna db" );

        if ($backendname eq 'Vmatch' ) {
            eval { $sbe->search({
                        query => 'GACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCA',
                        reverse_complement => 0,
                        upstream => 5,
                        showdesc => 100,
                  });  
            };
            ok( $EVAL_ERROR, 'Exception occured when upstream and showdesc' );
            eval { $sbe->search({
                        query => 'GACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCA',
                        reverse_complement => 0,
                        downstream => 5,
                        showdesc => 100,
                  });  
            };
            ok( $EVAL_ERROR, 'Exception occured when downstream and showdesc' );
            eval { $sbe->search({
                        query => 'GACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCA',
                        reverse_complement => 0,
                        upstream => 5,
                        downstream => 5,
                        showdesc => 100,
                  });  
            };
            ok( $EVAL_ERROR, 'Exception occured when up-&downstream and showdesc' );
            eval { $sbe->search({
                        query => 'GACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCA',
                        reverse_complement => 0,
                        showdesc => 100,
                  });  
            };
            ok( !$EVAL_ERROR, 'No exception occured without up/down' ) || diag $EVAL_ERROR;
            eval { $sbe->search({
                        query => 'GACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCA',
                        reverse_complement => 0,
                        complete => 1,
                        qspeedup => 2,
                  });  
            };
            ok( $EVAL_ERROR, 'No exception occured without up/down' );

        }    
        
################################################################################
        #  exit if $backendname eq "Hypa";
        # clean up
        CLEANUP:
        foreach my $file (@files) {
            $file = $sbe->is_word($file);
            unlink("t/data/$file");
        }
#        ok( rmdir("t/tmp"),
        #           "Can remove tmp directory (all temp files deleted)" );
        ok( rmdir("t/data"),
            "Can remove data directory (all data files deleted-just a test for the test)"
        );
    }    #SKIP
}

eval { Bio::Grep->new('UnknownBackend'); };
ok($EVAL_ERROR, 'Exception occured with unknown backend');

#ok( $tests > 0, " at least one backend found in path" );

# some helper functions
################################################################################
################################################################################
sub _get_ids_where_mm_smaller {
    my $mm = shift;
    my @results;
    foreach my $key ( keys %hits ) {
        push ( @results, $key ) if $hits{$key} <= $mm;
    }
    return sort @results;
}

sub compare_arrays {
    my ( $first, $second ) = @_;
    no warnings;    # silence spurious -w undef complaints
    return 0 unless @$first == @$second;
    for ( my $i = 0; $i < @$first; $i++ ) {
        return 0 if $first->[$i] ne $second->[$i];
    }
    return 1;
}

1;

# vim: ft=perl sw=4 ts=4 expandtab
