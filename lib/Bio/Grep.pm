package Bio::Grep;

use strict;
use warnings;

use Bio::Grep::Backend::Vmatch;
use Bio::Grep::Backend::Agrep;
use Bio::Grep::Backend::Hypa;
use Bio::Grep::Backend::GUUGle;
use Bio::Grep::Backend::RE;

use base 'Bio::Root::Root';

use version; our $VERSION = qv('0.8.2');

use Class::MethodMaker [
   new      => 'new2',
   scalar   => [qw / backend/],
];


sub new {
   my $self = shift->new2;
   my $backendname  = shift;
   $backendname = 'Vmatch' unless defined $backendname;
   my %known_backends = (
       Agrep  => 1,
       Vmatch => 1,
       Hypa   => 1,
       GUUGle => 1,
       RE     => 1,
   );    

   if (!defined $known_backends{$backendname}) {
      $self->throw(-class => 'Bio::Root::BadParameter',
                   -text  => 'Unknown back-end.',
                   -value => $backendname . ' not supported.');
   }

   my $backendpackage = "Bio::Grep::Backend::$backendname";
   
   $self->backend( $backendpackage->new() );
   $self->backend->settings->tmppath( File::Spec->tmpdir() );
   $self;
}
1;# Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep - Perl extension for searching in Fasta files

=head1 VERSION

This document describes Bio::Grep version 0.8.2

=head1 SYNOPSIS

  use Bio::Grep;
  
  my $search_obj = Bio::Grep->new('Vmatch');	
  
  # $sbe is now a reference to the back-end
  my $sbe = $search_obj->backend;
 
  # define the location of the suffix arrays
  $sbe->settings->datapath('data');
  
  mkdir($sbe->settings->datapath);	
  
  # now generate a suffix array. you have to do this only once.
  $sbe->generate_database('t/Test.fasta', 'Description for the test Fastafile');
  
  # search in this suffix array
  $sbe->settings->database('Test.fasta');
  
  # search for the reverse complement and allow 2 mismatches
  $sbe->settings->query('UGAACAGAAAG');
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(2);

  # or you can use Fasta file with queries
  # $sbe->settings->query_file('Oligos.fasta');

  # $sbe->search();

  # Alternatively, you can specify the settings in the search call.
  # This also resets everything except the paths and the database
  # (because it is likely that they don't change when search is called
  # multiple times)

  $sbe->search( { query  =>  'UGAACAGAAAG',
                  reverse_complement => 1,
                  mismatches         => 2,
                 });  
  
  my @ids;

  # output some informations! 
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->alignment_string() . "\n\n";
     push @ids, $res->sequence_id;
  }
  
  # get the gene sequences of all matches as Bio::SeqIO object.
  # (to generate a Fasta file for example)
  my $seqio = $sbe->get_sequences(\@ids);

=head1 DESCRIPTION

B<Bio-Grep> is a collection of Perl modules for searching in 
Fasta files. It supports different back-ends, most importantly some (enhanced) suffix
array implementations. Currently, there is no suffix array tool that works in
all scenarios (for example whole genome, protein and RNA data). Bio::Grep
provides a common API to the most popular tools. This way, you can easily switch or combine
tools.

=head1 METHODS

=over 

=item C<new($backend)>

This function constructs a C<Bio::Grep> object. Available external back-ends
are C<Vmatch>, C<Agrep>, C<GUUGle> and C<Hypa>. Perl regular
expressions are available in the C<RE> back-end. C<Vmatch> is default.

Sets temporary path to C<File::Spec-E<gt>tmpdir();>


  my $search_obj = Bio::Grep->new('Agrep');	


=item C<backend()>

Get/set the back-end. This is a object that uses L<Bio::Grep::Backend::BackendI> 
as base class. See L<Bio::Grep::Backend::BackendI>, L<Bio::Grep::Backend::Vmatch>,
L<Bio::Grep::Backend::Agrep>, L<Bio::Grep::Backend::GUUGle>,
L<Bio::Grep::Backend::RE> and L<Bio::Grep::Backend::Hypa>

=back

=head1 FEATURES

=over

=item 

Bio::Grep supports most of the features of the back-ends. If you need a 
particular feature that is not supported, write a feature request. In general it 
should be easy to integrate. For a complete list of supported features, see
L<Bio::Grep::SearchSettings>, for an overview see 
L<"FEATURE COMPARISON">.

=item 

This module should be suitable for large datasets. The back-end output is piped
to a temporary file and the parser only stores the current hit in memory.

=item

Bio::Grep has a nice interface for search result filters. See L<"FILTERS">.
This module also allows you to retrieve up- and downstream regions. Together
with filters, this makes Bio::Grep an ideal framework for seed and extend
algorithms. 

=item 

Bio::Grep was in particular designed for web services and therefore
checks the settings carefully before calling back-ends. See L<"SECURITY">.

=back

=head1 QUICK START

This is only a short overview of the functionality of this module.
You should also read L<Bio::Grep::Backend::BackendI> and the documentation of
the back-end you want to use (e.g. L<Bio::Grep::Backend::Vmatch>).

L<Bio::Grep::Cookbook> is a (not yet comprehensive) collection of recipes for
common problems.

=head2 GENERATE DATABASES 

As a first step, you have to generate a Bio::Grep database out of your Fasta
file in which you want to search. A Bio::Grep database consists of a couple of
files and allows you to retrieve informations about the database as well
as to perform queries as fast and memory efficient as possible. You have to do
this only once for every Fasta file.

For example:

  my $sbe = Bio::Grep->new('Vmatch')->backend;	
  $sbe->settings->datapath('data');
  $sbe->generate_database('../t/Test.fasta', 'Description for the test Fastafile');

Now, in a second script:

  my $sbe = Bio::Grep->new('Vmatch')->backend;	
  $sbe->settings->datapath('data');

  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
  

Alternatively, you can use L<bgrep> which is part of this distribution:

  bgrep --backend Vmatch --database TAIR6_cdna_20060907 --datapath data --createdb

  
=head2 SEARCH SETTINGS

All search settings are stored in the L<Bio::Grep::SearchSettings>
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

See L<Bio::Grep::SearchResult> for all available informations.


=head1 BGREP

This distribution comes with a sample script called B<bgrep>. See L<bgrep> for
details.

=head1 WHICH BACKEND?


We support these external back-ends:

=over

=item 

Vmatch (L<http://vmatch.de/>) 
	
=item 

Agrep (L<http://www.tgries.de/agrep/>)

=item 

GUUGle (L<http://bibiserv.techfak.uni-bielefeld.de/guugle/>)

=item 

Hypa (L<http://bibiserv.techfak.uni-bielefeld.de/HyPa/>) 

=back

=head2 FEATURE COMPARISON

=begin html

<table><tr><th>Feature</th><th>Agrep</th><th>GUUGle</th><th>HyPa</th><th>RE</th><th>Vmatch</th></tr><tr><td>Suffix Arrays/Trees</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Sliding Window</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
</tr>
<tr><td>Persistent Index<sup>1</sup></td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Mismatches</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Edit Distance</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Insertions</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
</tr>
<tr><td>Deletions</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
</tr>
<tr><td>Multiple Queries<sup>2</sup></td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>GU<sup>3</sup></td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
</tr>
<tr><td>DNA/RNA</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Protein</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Direct and Revcom</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Reverse Complement</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Upstream/Downstream</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Filters</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Query Length<sup>4</sup></td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
</tr>
<tr><td>Regular Expressions<sup>5</sup></td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
<td style="font-weight: bold;text-align: center;background-color: #00ff00;">yes</td>
<td style="text-align:center;background-color: #ffe0e0;">no</td>
</tr>
</table><br/><div style="font-size: smaller"><hr width="300" align="left"><sup>1</sup>Needs precalculation and (much) more memory but queries are in general faster<br/><sup>2</sup>With query_file<br/><sup>3</sup>HyPa also allows that GU counts only as 0.5 mismatches<br/><sup>4</sup>Matches if a substring of the query of size n or larger matches<br/><sup>5</sup>Agrep soon</div>

=end html

=begin man

     Features                 || Agrep  | GUUGle |  HyPa  |   RE   | Vmatch 
     Suffix Arrays/Trees      ||   no   |  yes   |  yes   |   no   |  yes   
     Sliding Window           ||  yes   |   no   |   no   |  yes   |   no   
     Persistent Index 1       ||   no   |   no   |  yes   |   no   |  yes   
     Mismatches               ||  yes   |   no   |  yes   |   no   |  yes   
     Edit Distance            ||  yes   |   no   |   no   |   no   |  yes   
     Insertions               ||   no   |   no   |  yes   |   no   |   no   
     Deletions                ||   no   |   no   |  yes   |   no   |   no   
     Multiple Queries 2       ||   no   |  yes   |   no   |   no   |  yes   
     GU 3                     ||   no   |  yes   |  yes   |   no   |   no   
     DNA/RNA                  ||  yes   |  yes   |  yes   |  yes   |  yes   
     Protein                  ||  yes   |   no   |  yes   |  yes   |  yes   
     Direct and Revcom        ||   no   |  yes   |   no   |  yes   |  yes   
     Reverse Complement       ||  yes   |  yes   |  yes   |  yes   |  yes   
     Upstream/Downstream      ||   no   |  yes   |  yes   |  yes   |  yes   
     Filters                  ||   no   |  yes   |  yes   |  yes   |  yes   
     Query Length 4           ||   no   |  yes   |   no   |   no   |  yes   
     Regular Expressions 5    ||   no   |   no   |   no   |  yes   |   no   

--
 1 Needs precalculation and (much) more memory but queries are in general faster
 2 With query_file
 3 HyPa also allows that GU counts only as 0.5 mismatches
 4 Matches if a substring of the query of size n or larger matches
 5 Agrep soon

=end man

Vmatch is fast but needs a lot of memory. Agrep is the best choice if you allow many 
mismatches in short sequences, if you want to search in Fasta files 
with relatively short sequences (e.g transcript databases) and if you are 
only interested in which sequences the approximate match was found. 
Its performance is in this case 
amazing. If you want the exact positions of a match in the sequence, choose vmatch. If you want 
nice alignments, choose vmatch too (EMBOSS can automatically align the 
sequence and the query in the agrep back-end, but then vmatch is faster). 
Filters require exact positions, so you can't use them with agrep. 
This may change in future version or not.

GUUGle may be the best choice if you have RNA queries (counts GU as no mismatch)
and if you are interested in only exact matches. Another
solution here would be to use Vmatch and write a filter (see next section) that
only allows GU mismatches. Of course, this is only an alternative if you
can limit (C<$sbe-E<gt>settings-E<gt>mismatches()>) the maximal number of GU
mismatches. Vmatch with its precalculated suffix arrays is really fast, so you
should consider this option.


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
      ... # here some local variables, see perldoc Class::MethodMaker
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
      ...
      $self->message('passed');
      return 1;
   }   
   
   sub reset {
      my $self = shift;
      # if you need local variables, you can clean up here
   }

   1;# Magic true value required at end of module

To apply your filter:

   ...

   my $filter = MyFilter->new();

   $sbe->settings->filters( ( $filter ) );
   $sbe->search();

See L<Bio::Grep::Filter::FilterI>.

=head1 ERROR HANDLING

Bio::Grep throws Bioperl exceptions when errors occur. You can use the 
module L<Error> to catch these exceptions:

   use Error qw(:try);
  
   ...
  
   try {
       $sbe->search();
   } catch Bio::Root::SystemException with {
       my $E = shift;
       print STDERR 'Back-end call failed: ' .     
       $E->{'-text'} . ' (' .  $E->{'-line'} . ")\n";
       exit(1);    
   } catch Bio::Root::BadParameter with {
       my $E = shift;
       print STDERR 'Wrong Settings: ' .     
       $E->{'-text'} . ' (' .  $E->{'-line'} . ")\n";
       exit(1);    
   } otherwise {        
       my $E = shift;
       print STDERR "An unexpected exception occurred: \n$E";
       exit(1);  
   };

Bio::Grep throws a C<SystemException> when a system() call failed,
C<BadParameter> whenever Bio::Grep recognizes some problems in the settings.
Be aware that Bio::Grep does not find all of these problems. In such a case,
the back-end call will fail and Bio::Grep will throw a C<SystemException>.

Whenever it is not possible to open, copy, close, delete or
write a file, croak() (L<Carp>) is called.

See L<Bio::Root::Exception>, L<Carp>.

=head1 SECURITY

The use of Bio::Grep (in Web Services for example) should be quite secure. All
test run in taint mode. Bio::Grep checks the settings before it generates the string
for the system() call. Bio::Grep uses L<File::Temp> for all temporary files.

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported. 

There is not yet a nice interface for searching for multiple queries. However,
Vmatch and GUUGle support this feature. So you can generate a Fasta query file
with L<Bio::SeqIO> and then set C<$sbe-E<gt>settings-E<gt>query_file()>. To
find out, to which query a match belongs, you have to check C<$res-E<gt>query>.

It is likely that C<$sbe-E<gt>settings-E<gt>query> is renamed to C<queries()>.

Please report any bugs or feature requests to
C<bug-bio-grep@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>. 


=head1 SEE ALSO

L<Bio::Grep::Cookbook>
L<Bio::Grep::Backend::BackendI>
L<Bio::Grep::Backend::Vmatch>
L<Bio::Grep::Backend::Agrep>
L<Bio::Grep::Backend::Hypa>
L<Bio::Grep::Backend::GUUGle>


=head2 PUBLICATIONS

GUUGle: L<http://bioinformatics.oxfordjournals.org/cgi/content/full/22/6/762>

HyPa:   L<http://nar.oxfordjournals.org/cgi/content/full/29/1/196>


=head1 AUTHOR

Markus Riester, E<lt>mriester@gmx.deE<gt>


=head1 LICENCE AND COPYRIGHT

Based on Weigel::Search v0.13

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
