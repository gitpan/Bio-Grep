#!perl -T

use Test::More;


BEGIN {
use lib 't';
use BioGrepTest;
my %prereq = BioGrepTest::check_prereq();
if (!$prereq{bioperl}) {
    plan skip_all => 'Bioperl not found';
}
elsif (!$prereq{bioperl_run}) {
    plan skip_all => 'Bioperl-run not found';
}
else {
    plan tests => 11;
}

use_ok( 'Bio::Grep' );
use_ok( 'Bio::Grep::Root' );
use_ok( 'Bio::Grep::Backends::Agrep' );
use_ok( 'Bio::Grep::Backends::BackendI' );
use_ok( 'Bio::Grep::Backends::GUUGle' );
use_ok( 'Bio::Grep::Backends::Hypa' );
use_ok( 'Bio::Grep::Backends::Vmatch' );
use_ok( 'Bio::Grep::Container::SearchResult' );
use_ok( 'Bio::Grep::Container::SearchSettings' );
use_ok( 'Bio::Grep::Filter::FilterI' );
use_ok( 'Bio::Grep::Filter::FilterRemoveDuplicates' );
}

diag( "Testing Bio::Grep $Bio::Grep::VERSION" );

