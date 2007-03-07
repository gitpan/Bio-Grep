package Bio::Grep;

use strict;
use warnings;

use Bio::Grep::Backends::Vmatch;
use Bio::Grep::Backends::Agrep;
use Bio::Grep::Backends::Hypa;
use Bio::Grep::Backends::GUUGle;

use base 'Bio::Root::Root';

use version; our $VERSION = qv('0.0.3');

use Class::MethodMaker [
   new      => 'new2',
   scalar   => [qw / backend/],
];


sub new {
   my $self = shift->new2;
   my $backendname  = shift;
   $backendname = 'Vmatch' unless defined $backendname;

   unless ($backendname eq 'Agrep' || $backendname eq 'Vmatch' || 
           $backendname eq 'Hypa'  || $backendname eq 'GUUGle') {
      $self->throw(-class => 'Bio::Root::BadParameter',
                   -text  => 'Unknown back-end.',
                   -value => $backendname . ' not supported.');
   }
   my $backendpackage = "Bio::Grep::Backends::$backendname";
   
   $self->backend( $backendpackage->new() );
   $self->backend->settings->tmppath( File::Spec->tmpdir() );
   $self;
}
1;# Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep - Perl extension for searching in Fasta files

=head1 VERSION

This document describes Bio::Grep version 0.0.3


=head1 SYNOPSIS

  use Bio::Grep;
  
  my $search_obj = Bio::Grep->new('Vmatch');	
  
  # $sbe is now a reference to the back-end
  # perldoc Bio::Grep::Backends::BackendI	
  # perldoc Bio::Grep::Backends::Vmatch
  my $sbe = $search_obj->backend;
 
  # perldoc Bio::Grep::Container::SearchSettings  
  $sbe->settings->datapath('data');
  
  mkdir($sbe->settings->datapath);	
  
  # generate a suffix array. you have to do this only once.
  $sbe->generate_database_out_of_fastafile('t/Test.fasta', 'Description for the test Fastafile');
  
  $sbe->settings->database('Test.fasta');
  
  # search for the reverse complement and allow 2 mismatches
  $sbe->settings->query('UGAACAGAAAG');
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(2);
  
  # $sbe->search();

  # Alternatively, you can specify the settings in the search call.
  # This also resets everything except the paths and the database
  # (because it is likely that they don't change when search is called
  # multiple times)

  $sbe->search( { query => 'UGAACAGAAAG',
                  reverse_complement => 1,
                  mismatches         => 2,
                 });  
  
  # output some informations! 
  # perldoc Bio::Grep::Container::SearchResult
  
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->alignment_string() . "\n\n";
  }


=head1 DESCRIPTION

Bio-Grep is a collection of Perl modules for searching in 
FASTA-files. It is programmed in a modular way. There are different 
back-ends available. You can filter search results.

=head1 METHODS

=over 

=item C<new($backend)>

This function constructs a C<Bio::Grep> object. Available back-ends
are Vmatch, Agrep, GUUGle and Hypa. Vmatch is default.

Sets temporary path to File::Spec->tmpdir();


  my $search_obj = Bio::Grep->new('Agrep');	


=item C<backend()>

Get/set the back-end. This is a object that uses L<Bio::Grep::Backend::BackendI> 
as base class. See L<Bio::Grep::Backends::BackendI>, L<Bio::Grep::Backends::Vmatch>,
L<Bio::Grep::Backends::Agrep>, L<Bio::Grep::Backends::GUUGle> and L<Bio::Grep::Backends::Hypa>

=back


=head1 QUICK START

This is only a short overview of the functionality of this module.
You should also read L<Bio::Grep::Backends::BackendI> and the documentation of
the back-end you want to use (e.g. L<Bio::Grep::Backends::Vmatch>).

=head2 GENERATE DATABASES 

As a first step, you have to generate a Bio::Grep database out of your FASTA
file in which you want to search. A Bio::Grep database consists of a couple of
files and allows you to retrieve informations about the database as well
as to perform queries as fast and memory efficient as possible. You have to do
this only once for every FASTA file.

For example:

  my $sbe = Bio::Grep->new('Vmatch')->backend;	
  $sbe->settings->datapath('data');
  $sbe->generate_database_out_of_fastafile('../t/Test.fasta', 'Description for the test Fastafile');

Now, in a second script:

  my $sbe = Bio::Grep->new('Vmatch')->backend;	
  $sbe->settings->datapath('data');

  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
  

=head2 SEARCH SETTINGS

All search settings are stored in the L<Bio::Grep::Container::SearchSettings>
object of the back-end: 

  $sbe->settings

To set an option, call

  $sbe->settings->optionname(value)

For example
    
  $sbe->settings->datapath('data');
  # take first available database 
  $sbe->settings->database($local_dbs[0]);

See the documentation of your back-end for available options. 

=head2 SEARCH

To start the back-end with the specified settings, simply call

  $sbe->search();

This method also accepts an hash reference with settings. In this case, all
previous defined options except all paths and the database are set to their
default values.

  $sbe->search({ mismatches => 2, 
                 reverse_complement => 0, 
                 query => $query });

=head2 ANALYZE SEARCH RESULTS

Use such a Bioperl like while loop to analyze the search results.

  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->alignment_string() . "\n\n";
  }

See L<Bio::Grep::Container::SearchResult> for all available informations.


=head1 BGREP

This distribution comes with a sample script called bgrep. See L<bgrep> for
details.

=head1 WHICH BACKEND?


We support this back-ends:

=over

=item Vmatch (L<http://vmatch.de/>) 
	
=item Agrep (L<http://www.tgries.de/agrep/>)

=item GUUGle (L<http://bibiserv.techfak.uni-bielefeld.de/guugle/>)

=item Hypa (L<http://bibiserv.techfak.uni-bielefeld.de/HyPa/>) 

=back

Vmatch is fast but needs a lot of memory. Agrep is the best choice if you allow many 
mismatches in short sequences, if you want to search in Fasta files 
with relatively short sequences (e.g transcript databases) and if you are 
only interested in which sequences the approximate match was found. 
Its performance is in this case 
amazing. If you want the exact positions, choose vmatch. If you want 
nice alignments, choose vmatch too (EMBOSS can automatically align the 
sequence and the query in the agrep back-end, but then vmatch is faster). 
Filters require exact positions, so you can't use them with agrep. 
This may change in future version or not.
GUUGle is the best choice if you have RNA queries (counts GU as no mismatch)
and if you are interested in only exact matches. 


=head1 FILTERS

Filtering search results is a common task. For that, Bio::Grep provides an
filter interface, L<Bio::Grep::Filter::FilterI>. Writing filters is
straightforward: 

   
   package MyFilter;
   
   use strict;
   use warnings;
   
   use Bio::Grep::Filter::FilterI;
   
   use base 'Bio::Grep::Filter::FilterI';
   
   use Class::MethodMaker
    [ new => [ qw / new2 / ],
       # here some local variables                 
    ];
   
   sub new {
      my $self = shift->new2;
      $self->delete(1); # a filter that actually filters, not only adds
                        # remarks to $self->search_result->remark

      $self->supports_alphabet( dna => 1, protein => 1);
      $self;
   }
   
   sub filter {
      my $self = shift;
      # code that examines $self->search_result
      # and returns 0 (not passed) or 1 (passed)
      return 1;
   }   
   
   sub reset {
      my $self = shift;
      # if you need local variables, can can clean it up here
   }

   1;# Magic true value required at end of module

To apply your filter:

   ...

   my $filter = MyFilter->new();

   $sbe->settings->filters( ( $filter ) );
   $sbe->search();

See L<Bio::Grep::Filter::FilterI>.

=head1 SECURITY

The use of Bio::Grep (in Web Services for example) should be quite secure. All
test run in taint mode. We check the settings before we generate the string
for the C<system()> call. We use L<File::Temp> for all temporary files.

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported. 

There is not yet a nice interface for searching for multiple queries. However,
Vmatch and GUUGle support this feature. So you can generate a FASTA query file
with L<Bio::SeqIO> and then set C<$sbe-E<gt>settings-E<gt>query_file()>. To
find out, to which query a match belongs, you have to check C<$res-E<gt>query>.

It is likely that C<$sbe-E<gt>settings-E<gt>query> is renamed to C<queries()>.

Please report any bugs or feature requests to
C<bug-bio-grep@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>. 


=head1 SEE ALSO

L<Bio::Grep::Backends::BackendI>
L<Bio::Grep::Backends::Vmatch>
L<Bio::Grep::Backends::Agrep>
L<Bio::Grep::Backends::Hypa>
L<Bio::Grep::Backends::GUUGle>


=head2 PUBLICATIONS

GUUGle: L<http://bioinformatics.oxfordjournals.org/cgi/content/full/22/6/762>

HyPa:   L<http://nar.oxfordjournals.org/cgi/content/full/29/1/196>


=head1 AUTHOR

Markus Riester, E<lt>mriester@gmx.deE<gt>


=head1 LICENCE AND COPYRIGHT

Based on Weigel::Seach v0.13

Copyright (C) 2005-2006 by Max Planck Institute for Developmental Biology, 
Tuebingen.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

# vim: ft=perl sw=4 ts=4 expandtab
