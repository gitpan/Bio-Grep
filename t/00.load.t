#!perl -T
use Test::More tests => 12;

use lib 't';

BEGIN {
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
use_ok( 'BioGrepTest' );
}

diag( "Testing Bio::Grep $Bio::Grep::VERSION" );

