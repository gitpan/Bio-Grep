#!perl -T
################################################################################
# some filter tests
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

use TestFilter;

my %tests = ( Agrep => 2, Vmatch => 14, GUUGle => 14, RE => 14 );
my $number_tests = 1;

foreach (keys %tests) {
    $number_tests += $tests{$_};
}
plan tests => $number_tests;

################################################################################

use Bio::Grep::Filter::FilterRemoveDuplicates;

use Bio::Grep;
use Bio::Perl;

# the number of files we assume in the data directory after
# database generation. NOT TESTED HERE
my %backend_filecnt = ( Agrep => 4, Vmatch => 15, GUUGle => 2  );

my $tests = 0;
mkdir("t/tmp");
mkdir("t/data");

foreach my $backendname ( sort keys %tests ) {
SKIP: {
        if ( $backendname ne 'RE' && BioGrepTest::find_binary_in_path( lc($backendname) ) eq '' ) {
            skip "$backendname not found in path", $tests{$backendname};
        }
        else {
            $tests++;
        }
        # diag("\n*** Testing $backendname ***");
        BioGrepTest::set_path( ( map { lc($_) } keys %backend_filecnt ),
            'RNAcofold' );
        my $gum = 1;
        $gum = 0 if $backendname eq 'GUUGle';
        my $sbe = Bio::Grep->new($backendname);

 # define own tmppath, so that we can check if all temporary files are deleted
        $sbe->settings->tmppath('t/tmp');
        $sbe->settings->datapath('t/data');
        BioGrepTest::delete_files;
        $sbe->generate_database( { file => 't/Test2.fasta',
            description => 'Description for Test2.fasta' });
        $sbe->settings->database('Test2.fasta');
        my $sequence = 'tgacagaagagagtgagcac';
        my $filter = Bio::Grep::Filter::FilterRemoveDuplicates->new();
        if ($backendname ne 'Agrep') {
            $sbe->search({
                    query   => 'ACCTAAAGTCACTG',
                    filters => [ $filter ],
                    reverse_complement => 0,
                    gumismatches => $gum,
                });
            my @ids;
            while (my $res = $sbe->next_res()) {
                my $rs_id  = $res->sequence->id;
                if ($rs_id =~ m{ \.\d \z }xms) {
                    ( $rs_id ) = $rs_id =~ m{ \A (.*?) \.\d \z }xms;
                }    
                push @ids, $rs_id;
            }
            is_deeply(\@ids, [ 'At2g42200' ], 'ids correct') || diag join(',',@ids);
        }    
        for my $j ( 0 .. 1 ) {
            if ($backendname eq 'RE') {
                $sbe->settings->query('[CG]TGC[AT]CTCTCTCTTCT[CG]TCA');
                $sbe->settings->reverse_complement(0);
            }   
            else {
                $sbe->settings->query($sequence);
                $sbe->settings->reverse_complement(1);
            }    
            if (defined $sbe->features->{MISMATCHES}) {
                $sbe->settings->mismatches(3);
            }
            if ($backendname eq 'GUUGle') {
                $sbe->settings->gumismatches(0);
                $sbe->settings->query_length(13);
            }    
            else {
                $sbe->settings->query_length_reset();
                $sbe->settings->gumismatches_reset();
            }
            if ( defined $sbe->features->{FILTERS}) {
                if ($j == 1) {
                    $filter->delete(0);
                    $sbe->settings->filters( $filter );
                } else {
                    $sbe->settings->filters_reset;
                }    
            }    
            $sbe->search();
            my @ids = BioGrepTest::get_sorted_result_ids($sbe);
#            diag join ',', @ids; 
            if ($backendname eq 'GUUGle') {
                is_deeply( \@ids, 
           [qw(At2g42200.1 At2g42200.1 At2g42200.2 At2g42200.2 At5g43270.1 At5g43270.2)],
                    "6 results" );
            }
            else {
                is_deeply( \@ids, 
                    [qw(At2g42200.1 At2g42200.2 At5g43270.1 At5g43270.2)],
                    "4 results" );
            }    
        }
        next unless defined( $sbe->features->{FILTERS} );
        $filter->delete(1);
        ok($filter->supports_alphabet_exists('dna'), 
           'Filter supports DNA');
        ok($filter->supports_alphabet_exists('protein'),
           'Filter supports Protein');
        ok(!$filter->supports_alphabet_exists('rna'),
           'Filter does not support RNA');
        my @filters = ( $filter );        
        for my $j ( 0 .. 1 ) {
            if ($j == 1 && defined $sbe->features->{FILTERS}) {
                push @filters, TestFilter->new();
                $sbe->settings->sort('ga');
            }    
            $sbe->settings->filters(@filters);
            $sbe->settings->query($sequence);
            $sbe->settings->reverse_complement(1);
            if (defined $sbe->features->{MISMATCHES}) {
                $sbe->settings->mismatches(3);
            }
            if ($backendname eq 'RE') {
                $sbe->settings->query('[CG]TGC[AT]CTCTCTCTTCT[CG]TCA');
                $sbe->settings->reverse_complement(0);
            }   
            elsif ($backendname eq 'GUUGle') {
                $sbe->settings->gumismatches(0);
                $sbe->settings->query_length(13);
            }
            else {
                $sbe->settings->query_length_reset();
                $sbe->settings->gumismatches_reset();
            }
            $sbe->search();
            my @ids = ();
            my $start_dg = -100;
            while (my $res = $sbe->next_res  ) {
                if ($j == 1 && defined $sbe->features->{FILTERS}) {
                    cmp_ok($res->dG, '>=', $start_dg, 'sorting works');
                    cmp_ok($res->dG, '<=', 0, 'sorting works');
                    is($res->remark,'passed', 'remark ok')
                }
                $start_dg = $res->dG;
                #  print STDERR "IS: " . $res->sequence->id . "\n";
                push ( @ids, $res->sequence->id );
            }
            @ids = sort @ids;
            is( @ids, 2, "2 results" );
            $sbe->settings->sort_reset;
            $sbe->settings->filters_reset;
        }
    }    #skip
}

foreach my $file (<t/data/*>) {
    my ($filename) = $file =~ m{ ( \A t/data/ [\d\w.\-\_]+ \z ) }xms;
    unlink $filename if -e $filename;
}

#ok( $tests > 0, " at least one backend found in path" );
#ok( rmdir("t/tmp"), "Can remove tmp directory (all temp files deleted)" );
ok( rmdir("t/data"),
    "Can remove data directory (all data files deleted-just a test for the test)"
);

# vim: ft=perl sw=4 ts=4 expandtab
