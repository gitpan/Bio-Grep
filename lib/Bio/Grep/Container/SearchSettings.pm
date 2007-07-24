package Bio::Grep::Container::SearchSettings;

use strict;
use warnings;

use Data::Dumper;

use version; our $VERSION = qv('0.7.0');

use Class::MethodMaker [
   new    => 'new2',
   scalar => [
      qw / mismatches insertions deletions editdistance query query_length
      _real_query gumismatches upstream downstream maxhits no_alignments
      datapath database online tmppath execpath reverse_complement sort
      complete query_file showdesc qspeedup hxdrop exdrop/
   ],
   array => [ qw / filters / ],
];

sub new {
   my $self = shift->new2;
    $self->_init();
   # initialize standard settings
   $self;
}

sub _init {
    my ( $self ) = @_;
   $self->mismatches(0);
   $self->insertions(0);
   $self->deletions(0);
   $self->editdistance_reset;
   $self->gumismatches(1);
   $self->no_alignments(0);
   $self->upstream(0);
   $self->downstream(0);
   $self->reverse_complement(0);
   $self->query_reset;
   $self->query_file_reset;
   $self->online_reset;
   $self->sort_reset;
   $self->complete_reset;
   $self->query_length_reset;
   $self->maxhits_reset;
   $self->showdesc_reset;
   $self->qspeedup_reset;
   return;
}    

sub set {
    my ( $self, $arg_ref ) = @_;
    $self->_init();
    for my $key (keys %$arg_ref) {
        $self->$key($arg_ref->{$key});
    } 
    return;
}    

sub to_string {
    my $self = shift;
    return Data::Dumper->Dump( [ $self ] ); 
}
1;# Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep::Container::SearchSettings - Data structure for all search settings 

=head1 SYNOPSIS

 
 use Bio::Grep::Backends::Vmatch;

 # configure our search back-end, in this case VMATCH
 my $sbe = Bio::Grep::Backends::Vmatch->new();
 
 # the back-end automatically generates a SearchSettings
 # object. You can use this now to configure the
 # back-end:  

 # tell the back-end the location of the executables
 # (only required if they are not in path)
 $sbe->settings->execpath('/usr/local/virtual.distrib');
 $sbe->settings->tmppath($tmp_path);
 $sbe->settings->mismatches(4);
 $sbe->settings->datapath($path_to_data_directory);
 $sbe->settings->database("ATH1.cdna");
 
 # retrieve upstream and downstream regions
 # see perldoc Bio::Grep::Container::SearchResult
 # how to get them out of the search results
 $sbe->settings->upstream(3);
 $sbe->settings->downstream(3);
 
 # search for the reverse complement 
 $sbe->settings->query("UGAACAGAAAGCUCAUGAGCC");
 $sbe->settings->reverse_complement(1);

 # Alternatively, you can specify the options in the search call:
 $sbe->search({
    query => $query,
    reverse_complement => 1,
 });

=head1 DESCRIPTION

B<Bio::Grep::Container::SearchSettings> is the data structure for all search settings.
Not all back-ends will support every option.

=head1 METHODS

=over 2

=item C<new()>

This function constructs an Bio::Grep::Container::SearchSettings object. The
back-end adds an object of this module into C<$sbe-E<gt>settings>, so you
should never have to call this constructor directly. If you want to reset all
settings to their default values, call

  $sbe->settings->set({});

See set().

=item C<set($hash_ref)>

Sets all settings in the hash reference:

  set( { query => $query, reverse_complement => 1 } );

This function resets everything except the paths and the database to default
values.

=item C<query()>

Get/Set the query, a simple string. Maybe this will change in future
versions to allow multiple queries.

   $sbe->settings->query('tgacagaagagagtgagcac');

=item C<query_file()>

Get/set query_file. The back-ends can create a query file based on the search
settings ($sbe->query()) or you can define one with $sbe->query_file.

    $sbe->settings->query_file('oligos.fasta');

Note that all settings that affect the creation of a query file (like reverse_
complement) may be ignored (depends on your back-end).

Currently only available in the Vmatch and GUUGle back-end.

=item C<reverse_complement()>

Get/set reverse_complement.

   #  search for the reverse complement
   $sbe->settings->reverse_complement(1)
   
   #  don't search for the reverse complement (default)
   $sbe->settings->reverse_complement(0)

=item C<mismatches()>

Get/Set allowed mismatches

   $sbe->settings->mismatches(5)

Not available in the GUUGle backend.

=item C<editdistance()>

Get/Set allowed edit distance. 

Only available in the Vmatch and Agrep back-end.

=item C<upstream()>

Get/set upstream. This is the number of bases upstream the match.
 
   $sbe->settings->upstream( 10 );

Not available in the Agrep backend.

=item C<downstream()>

Get/set downstream. This is the number of bases downstream the match.
   
   $sbe->settings->downstream( 10 );

Not available in the Agrep backend.

=item C<filters()>

Get/set the filters. This is an array of modules based on L<Bio::Grep::Filter::FilterI>.

   # display only possible targets of the miRNA query   
   my $filter1 = Bio::Grep::Filter::MIRNAFilter->new();
   
   # and display every gene id only once
   my $filter2 = Bio::Grep::Filter::FilterRemoveDuplicates->new();

   $sbe->settings->filters( ( $filter1, $filter2 ) );

Not available in the Agrep backend.

=item C<datapath()>

Get/set the datapath. This is the path to the databases for the back-end.

=item C<database()>

Get/set the database. This is the name of the database in $self->datapath

=item C<tmppath()>

Get/set the tmppath. This is a path were the back-end can store temporary files.

=item C<execpath()>

Get/set the execpath. This is a path to the back-end executable. If you don't
set it, make sure that the required executables are in path.

=item C<no_alignments()>

Get/set the no_alignments. Some back-ends like Agrep don't output alignments. 
The EMBOSS Smith-Waterman implementation will automatically generate alignments
for a search result. For perfomance reasons, you can turn turn this feature off:

  # if things need to be fast, turn alignments off
  $sbe->settings->no_alignments(1);

Not available (meaning not necessary) in the Vmatch and GUUGle back-end.

=item C<online()>

Get/set online. Vmatch back-end allows searching without using the index. When
you allow many mismatches, this could be faster. 

Only available in the Vmatch and HyPa back-end.

=item C<sort()>

Get set sort mode. You can get an hash with all available sort modes (including
a description from the back-end:

  my %sort_modes = $sbe->available_sort_modes();

=item C<maxhits()>

Get/set maxhits. Tells the back-end that it should output only the best n hits.

Only available in the Vmatch and GUUGle back-end.

=item C<gumismatches()>

Get/set gumismatches. Tells the back-end how it should count GU mismatches.
 Valid values in HyPa are 0, 0.5 and 1, in GUUGle only 0.

Only available in the HyPa and GUUGle back-end.

=item C<query_length()>

Get/Set the query length. Initialized with length of the query string. If
this member variable is smaller then the query string, then the back-end will
search for all substrings of that size. 

Only available in the Vmatch and GUUGle back-end.

=item C<to_string()> 

This function returns a string representation of this object


=back

=head2 Vmatch only

=over

=item C<complete()>

Get/set complete. Specify that query sequences must match completely.

    $sbe->settings->complete(1);

=item C<hxdrop()>

Specifies the xdrop value for hamming distance extension.

=item C<exdrop()>

Specifies the xdrop value for edit distance extension.

=item C<showdesc()>

Get/Set showdesc. This vmatch command line option makes the Vmatch parser
fetch the $sbe->result->sequence data directly out of the Vmatch output
instead of calling C<vseqsubselect>. Because of that, it is much faster with many
search results. You can't use this option if you want to retrieve up- or
downstream regions or if you are interested in the Vmatch internal sequence
id.

    # get the first 10 characters of the sequence description
    # (sequence id + annotation)
    $sbe->settings->showdesc(10);
    
=item C<qspeedup()>

Get/Set qspeedup. Specify speedup level when matching queries
(0: fast, 2: faster; default is 2). Beware of time/space tradeoff
             
=back

=head2 HyPa only

=over

=item C<insertions()>

Get/Set allowed insertions. 

=item C<deletions()>

Get/Set allowed deletions. 

=back


NOTE: You can use the hash C<features> from the back-end to check if some
feature is available or not. See L<Bio::Grep::Backends::BackendI>
for details.

=head1 SEE ALSO

L<Bio::Grep::Filter::FilterI>
L<Bio::Grep::Container::SearchResult>
L<Bio::Grep::Backends::BackendI>

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
