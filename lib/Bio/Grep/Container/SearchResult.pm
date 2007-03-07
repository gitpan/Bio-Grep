package Bio::Grep::Container::SearchResult;

use strict;
use warnings;

use IO::String;
use Carp qw(carp);

use base 'Bio::Root::Root';

use version; our $VERSION = qv('0.2.1');

use Class::MethodMaker [
   new    => 'new2',
   scalar => [qw / sequence query begin end alignment sequence_id remark percent_identity evalue dG _real_query/],
];

sub new {
   my $self = shift->new2;

   #initizialize member variables

   #from params
   $self->sequence(shift);
   $self->begin(shift);
   $self->end(shift);
   $self->alignment(shift);
   $self->sequence_id(shift);
   $self->remark(shift);
   $self;
}

sub mark_subject_uppercase {
   my $self   = shift;
   my $result = $self->sequence->seq;
   return
     lc( substr( $result, 0, $self->begin ) )
     . uc( substr( $result, $self->begin, $self->end - $self->begin ) )
     . lc( substr( $result, $self->end ) );
}

sub subject {
   my $self = shift;

   return Bio::Seq->new(
      -id  => $self->sequence->id, 
      -seq => $self->sequence->subseq( $self->begin + 1, $self->end ) );
}

sub alignment_string {
   my $self   = shift;
   my $result = "";
   my $str    = IO::String->new();
   my $out    = Bio::AlignIO->new( -format => 'clustalw', -fh => $str );
   unless (defined ($self->alignment)) {
      $self->warn("No alignment calculated.");
      return "";
   }   
   $out->write_aln( $self->alignment );
   $out->close();
   $str = IO::String->new( ${ $str->string_ref } );
   while ( my $l = <$str> ) {
      $result .= $l unless ( $l =~ /CLUSTAL/ or $l =~ /^\s+$/ );
   }
   return $result;
}

sub alignment_to_csv {

   # gets an alignment and returns and string for a CSV-file
   my $self = shift;

   my $aln = $self->alignment;
   unless (defined ($aln)) {
      $self->warn("No alignment calculated.");
      return "";
   }   
   return $self->get_gene_id() . ','
     . $_->_commafy_string( $aln->get_seq_by_pos(1)->seq ) . ","
     . $self->sequence->id . "\n,"
     . $self->_commafy_string( $aln->match_line )
     . ",\nQuery" . ","
     . $_->_commafy_string( $aln->get_seq_by_pos(2)->seq )
     . ",dG = "
     . $self->dG . "\n";
}

sub get_gene_id {
   my $self = shift;
   my $desc = $self->sequence->id;
   $desc =~ s/^(.*?)( |\|)//;
   return $1;
}

sub _commafy_string {

   # param: string
   # output: s,t,r,i,n,g
   my ( $self, $s ) = @_;
   if (!defined $s) {
      $self->throw(-class => 'Bio::Root::BadParameter',
                   -text  => 'Missing Parameter',);
   }
   
   my $result = "";
   for ( my $i = 0 ; $i != length($s) ; ++$i ) {
      $result .= substr( $s, $i, 1 ) . ',';
   }
   chop($result);
   return $result;
}
1;# Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep::Container::SearchResult - Data structure for a back-end search hit

=head1 SYNOPSIS
  
  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res ) {
     # $res->sequence is a Bio::Seq object with down-/upstream regions
     # see Bio::Grep::Container::SearchSettings
     print $res->sequence->id . "\n";
     
     # $res->subject is a Bio::Seq object without down-/upstream regions 
     print $res->subject->seq . "\n";

     # print down-/upstream regions lower case, subject sequence uppercase
     print $res->mark_subject_uppercase() . "\n";
     
     # output alignment
     print $res->alignment_string() . "\n";

     # print coordinates: perldoc Bio::SimpleAlign, Bio::LocatableSeq
     $print $res->alignment->get_seq_by_pos(1)->start . "\n\n";
  }

=head1 DESCRIPTION

This is the data structure for one hit in the database.

=head1 METHODS

=over 

=item C<new(sequence, begin, end, alignment, sequence_id, remark)>;

This function constructs a SearchResult object. 


=item C<sequence()>

Get/set the sequence found in database. L<Bio::Seq> object.

=item C<subject()>

Get the sequence found in database as string without upstream and downstream
regions. Bio::Seq object.

=item C<query()>

Get the query as L<Bio::Seq> object. Useful for multiple queries.

=item C<alignment()>

Get/set the alignment of the match. See L<Bio::SimpleAlign> for details. There are powerful
modules for converting this module in many formats. See L<Bio::AlignIO> for details.

=item C<get_gene_id()>

This function tries to guess gene name out of C<sequence-E<gt>name>.


=item C<sequence_id()>

Getter/Setter for sequence ID in database. This is an internal ID of the
back-end, not any ID of some annotation in the sequence name. Use
get_gene_id() to guess this ID. The internal ID can be used in the
back-end function get_sequences (L<Bio::Grep::Backends::BackendI>).

=item C<begin()>

Get/set the position of the beginning of the subject in the sequence. This allows 
retrieving upstream regions from the back-end.

=item C<end()>

Get/set the position of the end of the subject in the sequence. This allows retrieving
downstream regions from the back-end.

=item C<dG()>

Get/set dG . See L<Weigel::RNA::HybridizationI> for details.

=item C<remark()>

Get/set some additional informations like filter results to this hit.

=item C<evalue()>

Get/set the evalue of this hit.

=item C<percent_identity()>

Get/set the identity in percent of this hit. 

=back

=head1 OBJECT FORMATTER

Some predefined methods for printing objects. Make your own methods if you
don't like them. 

=over 

=item C<mark_subject_uppercase()>

This function returns the sequence in a string. the substring from $self->begin to
$self->end will be in uppercase, the rest in lowercase

=item C<alignment_string()>

This function returns a string with the formated alignment. We use CLUSTALW
Format without many blank lines and CLUSTAL header. In some back-ends like
Agrep, this function will return an empty string if no_alignments is true
in the back-end search settings (L<Bio::Grep::Container::SearchSettings>).

=item C<alignment_to_csv()>

This function generates a string for a CSV file, where every base of the stored alignment gets its own column. Additional infos
like description and dG are displayed, too.

=back

=head1 SEE ALSO

L<Bio::SimpleAlign> 
L<Bio::LocatableSeq> 
L<Bio::AlignIO> 
L<Bio::Seq>
L<Bio::SeqIO>
L<Weigel::RNA::HybridizationI>
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
