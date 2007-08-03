package Bio::Grep::Backend::RE;

use strict;
use warnings;

use Fatal qw(open close);

use Bio::Grep::SearchResult;
use Bio::Grep::Backend::BackendI;

use base 'Bio::Grep::Backend::Agrep';

use version; our $VERSION = qv('0.8.4');

sub new {
    my $self = shift;
    $self = $self->SUPER::new;
    my %all_features = $self->features;

    $all_features{DIRECT_AND_REV_COM} = 1;
    $all_features{UPSTREAM} = 1;
    $all_features{DOWNSTREAM} = 1;
    $all_features{FILTERS} = 1;
    $all_features{SORT} = 1;
    delete $all_features{MISMATCHES};
    delete $all_features{DELETIONS};
    delete $all_features{INSERTIONS};
    delete $all_features{EDITDISTANCE};
    $self->features(%all_features);
    $self;
}

sub search {
    my ($self, $arg_ref) = @_;
    my $s    = $self->settings;
    $self->_check_search_settings($arg_ref);
    my $regex = $s->query; 
    if (!defined $regex) {
        $self->throw( -class => 'Bio::Root::BadParameter',
                      -text  => 'Query not defined.',
                    );
    }    

    if ($s->direct_and_rev_com || $s->reverse_complement) {
        $self->throw( 
             -class => 'Bio::Root::BadParameter',
             -text => 'Query does not look like a DNA/RNA sequence.', 
             -value => $regex,
         ) if $regex !~ m{\A [gactu]+ \z}xmsi;
        my $tmp = Bio::Seq->new(-seq => $regex); 
        if ($s->direct_and_rev_com) {
            $regex = $regex .'|'.  $tmp->revcom->seq;
        }
        else {
            $regex =  $tmp->revcom->seq;
        }    
    }
    
    # I need the prefix to calculate the position of
    # the regex match
    $self->{_regex} = qr{$regex}imsx;
    $self->{_puffer} = [];
    

    open my $FH, '<', $self->_cat_path_filename( $s->datapath, $s->database . '.dat' );
    $self->_output_fh($FH);

    $self->_load_mapping();
    $self->_prepare_results;
    return 1;
}


sub _parse_next_res {
    my $self    = shift;
    my $s       = $self->settings;
    my @puffer = @{ $self->{_puffer} };
    if (scalar @puffer > 0) {
        return shift @{ $self->{_puffer} };
    }    

    my %mapping = $self->_mapping;

    my $FH = $self->_output_fh;
    while (my $line = <$FH>) {
        chomp $line;
        my ( $pos, $complete_seq ) = $line =~ m{\A (\d+) : (.*) \z}xms;

        # store sequence for upstream/downstream
        my $seq = $complete_seq;

        my $lastpos;
        while ( $seq =~ /$self->{_regex}/g ) {
            my $subject_pos = length $`;
            my $subject_seq = $&;
            $lastpos = $subject_pos + length($subject_seq); 
            my $up_seq = '';
            my $down_seq = '';

            #warn "$subject_seq";
            if ($s->upstream > 0 || $s->downstream > 0) {
                my $upstream_pos = $subject_pos - $s->upstream;
                $upstream_pos = 0 if $upstream_pos < 0;
                $up_seq = substr $complete_seq, $upstream_pos, $subject_pos -
                    $upstream_pos;
                    #warn "$upstream_pos   $subject_pos  $up_seq";

                my $downstream_pos = $lastpos+$s->downstream;
                $downstream_pos = length $complete_seq if $downstream_pos >
                length $complete_seq;
                $down_seq = substr $complete_seq, $lastpos, $downstream_pos -
                    $lastpos;
                    # warn "$downstream_pos   $lastpos  $down_seq";
            }    
            my $upstream = length $up_seq;    
            my $fasta = Bio::Seq->new(
                -id  => $mapping{$pos},
                -seq => $up_seq .  $subject_seq . $down_seq,
            );
            my $tmp_aln = new Bio::SimpleAlign( -source => "Bio::Grep" );
            $tmp_aln->add_seq(
                Bio::LocatableSeq->new(
                    -id    => 'Subject',
                    -seq   => $subject_seq,
                    -start => $subject_pos+1,
                    -end   => $subject_pos+ length($subject_seq),
                )
            );
            my $res = $self->_filter_result(
                Bio::Grep::SearchResult->new( $fasta,
                $upstream, $upstream +length($subject_seq),
                $tmp_aln, $fasta->id, '' )
                );
            if ($res) {    
                $res->query(Bio::Seq->new(-id => 'Query'));    
                push @puffer, $res;
            }

        }    
        if (scalar @puffer > 0) {
            $self->{_puffer} = \@puffer;
            return shift @{ $self->{_puffer} };
        }
    }
    return;
}

sub available_sort_modes {
    my ($self) = @_;
    return (
        ga => 'ascending order of dG',
        gd => 'descending order of dG',
    );
}

1;    # Magic true value required at end of module
__END__


=head1 NAME

Bio::Grep::Backend::RE - Perl Regular Expression back-end  


=head1 SYNOPSIS

  use Bio::Grep;
  
  my $sbe = Bio::Grep->new('RE');
  
  $sbe->settings->datapath('data');
  
  # generate a database. you have to do this only once. 
  $sbe->generate_database('ATH1.cdna', 'AGI Transcripts (- introns, + UTRs)');
  
  # search on both strands  
  # retrieve up- and downstream regions of size 30
  
  $sbe->search({
    query   => 'GAGCCCTT',
    direct_and_rev_com => 1, 
    upstream           => 30,
    downstream         => 30,
    database           => 'ATH1.cdna',
  });
  
  my @internal_ids;
  
  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";
     push @internal_ids, $res->sequence_id;
  }
  
  # get the complete sequences as Bio::SeqIO object
  my $seq_io = $sbe->get_sequences(\@internal_ids);

  # sequences with at least 10 As
  $sbe->search({ query => '[A]{10,}' });
 
  # some SNPs
  $sbe->search({query => '[CG]TGC[AT]CTCTTCT[CG]TCA'});

=head1 DESCRIPTION

B<Bio::Grep::Backend::RE> searches for a query with a
Perl Regular Expression. 

Note 1: B<BETA RELEASE!> 

Note 2: C<reverse_complement> (and C<direct_and_rev_com> ) are supported, but are
only available for DNA/RNA queries, not for regular expressions.

Note 3: The i modifier is added to the regex. This means that the search is
case insensitive. 

Note 4: Be careful with RNA sequences: U is not the same as T in this back-end!

Note 5: L<Bio::Grep::Backend::RE> databases are compatible with
L<Bio::Grep::Backend::Agrep> databases.

=head1 METHODS

See L<Bio::Grep::Backend::BackendI> for inherited methods. 

=over 2

=item C<Bio::Grep::Backend::RE-E<gt>new()>

This function constructs an RE back-end object

   my $sbe = Bio::Grep::Backend::RE->new();

=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

   $sbe->sort('ga');

Available sortmodes in RE:

=over

            ga  : 'ascending order of dG'
            gd  : 'descending order of dG'

=back

Note that 'ga' and 'gd' require that search results have dG set. 
L<Bio::Grep::RNA> ships with filters for free energy calculation. Also note that
these two sort options require that we load all results in memory.

=back

=head1 DIAGNOSTICS

See L<Bio::Grep::Backend::BackendI> for other diagnostics. 

=over

=item C<Query does not look like a DNA/RNA sequence.>

Either C<reverse_complement> or C<direct_and_rev_com> is set and the query
does not match the regular expression C<m{\A [gactu]+ \z}xmsi>. C<Bio::Root::BadParameter>.

=item C<Query not defined.>

You forgot to define C<$sbe-E<gt>settings-E<gt>query>. C<Bio::Root::BadParameter>.

=back

=head1 SEE ALSO

L<Bio::Grep::Backend::BackendI>
L<Bio::Grep::SearchSettings>
L<Bio::SeqIO>
L<Bio::Index::Fasta>

=head1 AUTHOR

Markus Riester, E<lt>mriester@gmx.deE<gt>

=head1 LICENCE AND COPYRIGHT

Copyright (C) 2007  by M. Riester. All rights reserved. 

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
