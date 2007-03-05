#!perl -T
################################################################################
# some tests for helper functions
#
################################################################################

use Test::More tests => 31;
################################################################################

use strict;
use warnings;

use English qw( -no_match_vars );
use Cwd;
use Scalar::Util qw/tainted/;

use Bio::Grep;
use Bio::Perl;

use lib 't';
use BioGrepTest;

my @paths = ( '', '/', '/usr/local/bin' );

my $sbe = Bio::Grep->new()->backend;

my $result = Bio::Grep::Container::SearchResult->new();

# todo make this platform independent
is( $sbe->_cat_path_filename( $paths[0], 't.txt' ), 't.txt', 'concat path' );

is( $result->_commafy_string('string'), 's,t,r,i,n,g' );
is( $result->_commafy_string(''),       '' );
eval { $result->_commafy_string(); };
ok($EVAL_ERROR);

my $tainted_word    = 'bla' . substr( cwd, 0, 0 );
my $tainted_integer = '1' . substr( cwd,   0, 0 );
my $tainted_real    = '1.1' . substr( cwd, 0, 0 );

ok( tainted $tainted_word,    $tainted_word . ' tainted' );
ok( tainted $tainted_integer, $tainted_integer . ' tainted' );
ok( tainted $tainted_real,    $tainted_real . ' tainted' );

my $not_tainted_integer = $sbe->is_integer($tainted_integer);
ok( !tainted $not_tainted_integer, $not_tainted_integer . ' not tainted' );
my $not_tainted_word = $sbe->is_word($tainted_word);
ok( !tainted $not_tainted_word, $not_tainted_word . ' not tainted' );
my $not_tainted_real = $sbe->is_real($tainted_real);
ok( !tainted $not_tainted_real, $not_tainted_real . ' not tainted' );

is( $sbe->is_integer('1234'), 1234 );
eval { $sbe->is_integer('1234.5'); };
ok($EVAL_ERROR);
eval { $sbe->is_integer('10 && ls *'); };
ok($EVAL_ERROR);
is( $sbe->is_integer(undef), undef );

is( $sbe->is_real('1234'),     1234 );
is( $sbe->is_real('1234.'),    '1234.' );
is( $sbe->is_real('-1234.12'), '-1234.12' );
is( $sbe->is_real(undef),      undef );
eval { $sbe->is_real('1234.5.1'); };
ok($EVAL_ERROR);
eval { $sbe->is_real('1234 .5'); };
ok($EVAL_ERROR);
eval { $sbe->is_real('10 && ls *'); };
ok($EVAL_ERROR);

is( $sbe->is_word('1234'),           1234 );
is( $sbe->is_word('1234-valid.txt'), '1234-valid.txt' );
is( $sbe->is_word('1234-valid.txt_'), '1234-valid.txt_' );
eval { $sbe->is_word('valid && ls *'); };
ok($EVAL_ERROR);

$sbe=Bio::Grep->new('GUUGle')->backend;

ok($sbe->_rnas_match('agcua','agcua'), 'rna matching function');
ok(!$sbe->_rnas_match('agcuag','agcua'), 'rna matching function');
ok($sbe->_rnas_match('agcuag','aucuag'), 'rna matching function');
ok($sbe->_rnas_match('agcuau','aucuag'), 'rna matching function');
ok($sbe->_rnas_match('aucuau','agcuag'), 'rna matching function');
ok(!$sbe->_rnas_match('aucaau','agcuag'), 'rna matching function');

# vim: ft=perl sw=4 ts=4 expandtab
