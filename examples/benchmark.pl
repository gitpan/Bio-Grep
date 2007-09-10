#!/usr/bin/perl
use strict;
use warnings;

use Time::HiRes qw(gettimeofday tv_interval);

use Bio::Grep;
use Bio::Perl;
use Data::Dumper;
use English qw( -no_match_vars ) ; 
use Template;

my $template = Template->new();

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
my $filenameCDNA = 'TAIR7_cdna_20070425';
my $iterations = 20;
my $iterationsdb = 2;

DB:
for $b (sort keys %be) {
    my $sbe = $be{$b};
    $time = [gettimeofday];
    for my $i (1..$iterationsdb) {
        system("rm -rf data$b/");
        mkdir 'data' . $b;

        $sbe->generate_database({
                datapath      => 'data' . $b,
                file          => "examples/$filenameCDNA",
                prefix_length => 3,
            }); 
    }
    $results{"${b}_dbgen"} = sprintf("%.2f",
        (tv_interval($time)/$iterationsdb));
    warn "$b took " . $results{"${b}_dbgen"} . " seconds\n";
}    
#goto CREATETMP;

MM:
for $b (sort keys %be) {
    my $sbe = $be{$b};
    my $loop_counter = 0;
    $loop_counter = 1 if $b eq 'vmatch';
    for my $online ( 0 .. $loop_counter) {
    for my $mm (0..5) {
        next MM if !defined $sbe->features->{MISMATCHES} && $mm > 0;
        $time = [gettimeofday];
        for my $i (1..$iterations) {
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
            online             => $online,
            database           => $filenameCDNA,
            gumismatches       => $gu,
            %showdesc,
            }); };
            my @ids;
            while (my $res = $sbe->next_res) {
                push @ids, $res->sequence->id;
            }    
            warn scalar(@ids). " results.\n" if $VERBOSE;
        }   
        $results{"${b}_mm_${mm}_$online"} = sprintf("%.2f",
            tv_interval($time)/$iterations);
        warn "$b (mm $mm) took " . $results{"${b}_mm_${mm}_$online"} . " seconds\n";
    }     
    }
}    

CREATETMP:

$results{cpuinfo} = get_cpu();
$results{filenameCDNA} = $filenameCDNA;
$results{biogrepv} = $Bio::Grep::VERSION;
$results{iterations} = $iterations;
$results{iterationsdb} = $iterationsdb;
$template->process('examples/Benchmarks.tt', \%results, 'lib/Bio/Grep/Benchmarks.pod') || die
$template->error(), "\n";

print Dumper \%results;

sub get_cpu {
    use Sys::Info::CPU;
    return scalar Sys::Info::CPU->new->identify;    
}    
