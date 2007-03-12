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

my %tests = ( Agrep => 3, Vmatch => 9, Hypa => 9, GUUGle => 9 );
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
my %backend_filecnt = ( Agrep => 5, Vmatch => 16, Hypa => 31, GUUGle => 3  );

my $tests = 0;
mkdir("t/tmp");
mkdir("t/data");

foreach my $backendname ( sort keys %tests ) {
SKIP: {
        if ( BioGrepTest::find_binary_in_path( lc($backendname) ) eq '' ) {
            skip "$backendname not found in path", $tests{$backendname};
        }
        else {
            $tests++;
        }
        BioGrepTest::set_path( ( map { lc($_) } keys %backend_filecnt ),
            'RNAcofold' );
        my $search_obj = Bio::Grep->new($backendname);
        my $sbe        = $search_obj->backend;

 # define own tmppath, so that we can check if all temporary files are deleted
        $sbe->settings->tmppath('t/tmp');
        $sbe->settings->datapath('t/data');
        BioGrepTest::delete_files;

        $sbe->generate_database_out_of_fastafile( 't/Test2.fasta',
            'Description for Test2.fasta' );

        $sbe->settings->database('Test2.fasta');
        my $sequence = 'tgacagaagagagtgagcac';
        for my $j ( 0 .. 2 ) {
            $sbe->settings->query($sequence);
            $sbe->settings->reverse_complement(1);
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
            $sbe->search();
            my @ids = ();
            foreach my $res ( @{ $sbe->results } ) {

                #   print STDERR "IS: " . $res->sequence->id . "\n";
                push ( @ids, $res->sequence->id );
            }
            @ids = sort @ids;
            if ($backendname eq 'GUUGle') {
                is( @ids, 6, "6 results" );
            }
            else {
                is( @ids, 4, "4 results" );
            }    
        }
        next unless defined( $sbe->features->{FILTERS} );
        my $filter = Bio::Grep::Filter::FilterRemoveDuplicates->new();
        ok($filter->supports_alphabet_exists('dna'));
        ok($filter->supports_alphabet_exists('protein'));
        ok(!$filter->supports_alphabet_exists('rna'));
        
        for my $j ( 0 .. 2 ) {
            $sbe->settings->filters($filter);
            $sbe->settings->query($sequence);
            $sbe->settings->reverse_complement(1);
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
            $sbe->search();
            my @ids = ();
            while (my $res = $sbe->next_res  ) {

                #  print STDERR "IS: " . $res->sequence->id . "\n";
                push ( @ids, $res->sequence->id );
            }
            @ids = sort @ids;
            is( @ids, 2, "2 results" );
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
