#!/usr/bin/perl
use strict;
use warnings;

use Time::HiRes qw(gettimeofday tv_interval);

use Bio::Grep;
use Bio::Perl;
use Data::Dumper;
use English qw( -no_match_vars ) ; 
my %be = (
   agrep  => Bio::Grep->new('Agrep'),
   vmatch => Bio::Grep->new('Vmatch'),
   re     => Bio::Grep->new('RE'),
   guugle => Bio::Grep->new('GUUGle'),
);

my %results;

my $query = 'ugacagaagagagugagcac';
$query =~ tr{u}{t};

my $time;
my $VERBOSE=0;

DB:
for $b (sort keys %be) {
    my $sbe = $be{$b};
    mkdir 'data' . $b;
    $time = [gettimeofday];
    $sbe->generate_database({
            datapath      => 'data' . $b,
            file          => 'testbig.fasta',
            prefix_length => 3,
        }); 
    $results{$b}{dbgen} = tv_interval($time);
    warn "$b took " . $results{$b}{dbgen} . " seconds\n";
}    

MM:
for $b (sort keys %be) {
    my $sbe = $be{$b};
    for my $mm (0..5) {
        next MM if !defined $sbe->features->{MISMATCHES} && $mm > 0;
        $time = [gettimeofday];
        for my $i (1..10) {
            my %showdesc;
            %showdesc = ( showdesc => 100) if $b eq 'vmatch'; 
            my $gu = 1;
            $gu = 0 if $b eq 'guugle';
            eval { $sbe->search({
            query              => $query,
            mismatches         => $mm,
            reverse_complement => 1,
            datapath => 'data' . $b,
            no_alignments      => 1,
            database           => 'testbig.fasta',
            gumismatches       => $gu,
            %showdesc,
            }); };
            my @ids;
            while (my $res = $sbe->next_res) {
                push @ids, $res->sequence->id;
            }    
#           print join (',',@ids). "\n";
            warn scalar(@ids). " results.\n" if $VERBOSE;
        }   
        $results{$b}{"mm$mm"} = tv_interval($time)/10;
        warn "$b (mm $mm) took " . $results{$b}{"mm$mm"} . " seconds\n";
    }     
}    

print Dumper \%results;
