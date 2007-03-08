#!perl

use strict;
use warnings;
use English qw(-no_match_vars);
use Test::More;

eval {
    require Test::Perl::Critic;
    my @config = (); #Arguments for Perl::Critic->new() go here!
    Test::Perl::Critic->import( @config );
};

if( $EVAL_ERROR ) {
    plan( skip_all => 'Test::Perl::Critic required for PBP tests' );
}

Test::Perl::Critic::all_critic_ok();
