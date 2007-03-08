#!perl -T
use Test::More tests => 1;

use lib 't';

BEGIN {
use_ok( 'BioGrepTest' );
}
use BioGrepTest;

my %prereq = BioGrepTest::check_prereq();
my $bp =<<EOT
* Bioperl    : NOT FOUND                                                    *
EOT
;
my $bp_run =<<EOT
* Bioperl-run: NOT FOUND                                                    *
EOT
;
my $backend = <<EOT
* Backend    : No backend found in path. You should install the backends    *
*              before running these tests. This way you make sure that the  *
*              parsers work with your installed version of the backend.     *
EOT
;

$bp =~ s/NOT FOUND/FOUND    / if $prereq{bioperl};
$bp_run =~ s/NOT FOUND/FOUND    / if $prereq{bioperl_run};
$backend = '' if $prereq{backend};

diag("\n");
diag('*' x 77);
diag($bp);
diag($bp_run);
diag($backend) if $backend ne '';
diag('*' x 77);

1;
