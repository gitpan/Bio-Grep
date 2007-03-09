#!perl -T
################################################################################
# some tests for helper functions
#
################################################################################

################################################################################

use strict;
use warnings;
BEGIN{
    use Test::More;
    use lib 't';
    use BioGrepTest;
    my %prereq = BioGrepTest::check_prereq();
    if (!$prereq{bioperl}) {
        plan skip_all => 'Bioperl not found';
    }
    elsif (!$prereq{bioperl_run}) {
        plan skip_all => 'Bioperl-run not found';
    }
}

plan tests => 32;

use English qw( -no_match_vars );
use Cwd;
use Data::Dumper;
use Scalar::Util qw/tainted/;

use Bio::Grep;
use Bio::Perl;


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
ok($sbe->_rnas_match('uguggu','cgcgau'), 'rna matching function');
ok($sbe->_rnas_match('uguggu','ugcggu'), 'rna matching function');
ok($sbe->_rnas_match('uguggu','cguggu'), 'rna matching function');
ok(!$sbe->_rnas_match('uguggu','cgcguu'), 'rna matching function');

my $tmp = $sbe->settings->tmppath;
$sbe->settings->datapath('data');
$sbe->settings->database('Test.fasta');
$sbe->settings->reverse_complement(1);

my $settings_dump =<<EOT
\$VAR1 = bless( {                               
                 'datapath' => 'data',
                 'no_alignments' => 0,
                 'execpath' => '',
                 'database' => 'Test.fasta',
                 'deletions' => '0',
                 'upstream' => '0',
                 'insertions' => '0',
                 'reverse_complement' => 1,
                 'tmppath' => '$tmp',
                 'mismatches' => '',
                 'downstream' => '0',
                 'gumismatches' => 0
               }, 'Bio::Grep::Container::SearchSettings' );
EOT
;

is_deeply(d2h($sbe->settings->to_string), d2h($settings_dump), 'Settings dump ok');
sub d2h {
    my ( $dump ) = @_;
    my %h;
    while ( $dump =~ m{ ^ \s+ '(.*?)' .*? > \s (.*?) [,]* $ }xmsg ) {
        my ($v1, $v2) = ($1, $2);
        $v2 =~ s/\'//g;
        $v2 = '' if !$v2;
        chomp $v2;
        $h{$v1} = $v2;
    }
    return \%h;
}

# vim: ft=perl sw=4 ts=4 expandtab
