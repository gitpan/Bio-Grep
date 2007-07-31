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
    plan tests => 12;
}

use_ok( 'Bio::Grep' );
use_ok( 'Bio::Grep::Root' );
use_ok( 'Bio::Grep::Backend::Agrep' );
use_ok( 'Bio::Grep::Backend::BackendI' );
use_ok( 'Bio::Grep::Backend::GUUGle' );
use_ok( 'Bio::Grep::Backend::Hypa' );
use_ok( 'Bio::Grep::Backend::RE' );
use_ok( 'Bio::Grep::Backend::Vmatch' );
use_ok( 'Bio::Grep::SearchResult' );
use_ok( 'Bio::Grep::SearchSettings' );
use_ok( 'Bio::Grep::Filter::FilterI' );
use_ok( 'Bio::Grep::Filter::FilterRemoveDuplicates' );
}

diag( "Testing Bio::Grep $Bio::Grep::VERSION" );

