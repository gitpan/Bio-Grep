package Bio::Grep::Backend::Agrep;

use strict;
use warnings;

use Fatal qw(open close);

use Bio::Grep::SearchResult;
use Bio::Grep::Backend::BackendI;

use base 'Bio::Grep::Backend::BackendI';

use File::Basename;
use IO::String;

use version; our $VERSION = qv('0.8.3');

sub new {
    my $self = shift;
    $self = $self->SUPER::new;
    my %all_features = $self->features;

    # agrep does not supprt this features
    delete $all_features{NATIVE_ALIGNMENTS};
    delete $all_features{GUMISMATCHES};
    delete $all_features{FILTERS};
    delete $all_features{EVALUE};
    delete $all_features{PERCENT_IDENTITY};
    delete $all_features{UPSTREAM};
    delete $all_features{DOWNSTREAM};
    delete $all_features{ONLINE};
    delete $all_features{SORT};
    delete $all_features{MAXHITS};
    delete $all_features{COMPLETE};
    delete $all_features{QUERY_FILE};
    delete $all_features{QUERY_LENGTH};
    delete $all_features{SHOWDESC};
    delete $all_features{QSPEEDUP};
    delete $all_features{HXDROP};
    delete $all_features{EXDROP};
    delete $all_features{REVCOM_DEFAULT};
    delete $all_features{DIRECT_AND_REV_COM};
    delete $all_features{NATIVE_D_A_REV_COM};
    $self->features(%all_features);
    $self;
}

sub search {
    my ($self, $arg_ref) = @_;
    my $s    = $self->settings;
    $self->_check_search_settings($arg_ref);
    
    my $query = $self->_prepare_query();

    # now generate the command string
    my $mm = 0;
    if ($s->mismatches_isset && $s->mismatches > 0) {
        $mm = $s->mismatches;
    }    

    # make insertions and deletions to expensive
    my $fuzzy = "-$mm -D" . ( $mm + 1 ) . " -I" . ( $mm + 1 );

    if ( $s->editdistance_isset ) {
        $fuzzy = "-" . $s->editdistance;
    }


    my $command =
        $self->_cat_path_filename( $s->execpath, 'agrep' ) . ' -i ' . $fuzzy
        . ' '
        . $query . ' '
        . $self->_cat_path_filename( $s->datapath, $s->database . '.dat' );

    if ( $ENV{BIOGREPDEBUG} ) {
        warn $command . "\n";
    }

    my $cmd_ok = $self->_execute_command($command);

    # delete temporary files
    #unlink($tmp_query_file) if !$query_file;

    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  => "Agrep call failed. Command was:\n\t$command"
        )
        if !$cmd_ok;
    #$self->_output_parser($output);
    $self->_load_mapping();
    $self->_prepare_results;
    return 1;
}

sub get_databases {
    my $self = shift;
    return $self->_get_databases('.map');
}

sub generate_database {
    my ( $self, $file, $description ) = @_;
    my ( $filename ) = fileparse($file);

    $self->_copy_fasta_file_and_create_nfo( $file, $filename, $description );

    $filename = $self->_cat_path_filename($self->settings->datapath,
        $filename);

    open my $DATFILE, '>',  "$filename.dat"; 
    
    open my $MAPFILE, '>', "$filename.map";

    my $in = Bio::SeqIO->new( -file => $filename, -format => 'Fasta' );
    my $i = 0;
    
    while ( my $seq = $in->next_seq() ) {
        print $MAPFILE $i . "\t" . $seq->id
            . "\n";
        print $DATFILE $i . ':' . $seq->seq
            . "\n";
        $i++;
    }
    close $DATFILE;
    close $MAPFILE;
    $self->_create_index_and_alphabet_file( $filename );
    return $filename;    
}

sub _load_mapping {
    my ( $self ) = @_;
    my $s       = $self->settings;
    my $mapfile   = $s->datapath . '/' . $s->database . '.map';


    my %mapping = ();

    open my $MAPFILE, '<', $mapfile; 

    while ( my $line = <$MAPFILE> ) {
        chomp $line;
        my ( $pos, $id ) = split "\t", $line;
        $mapping{$pos} = $id;
    }
    close $MAPFILE;

    $self->_mapping( %mapping );
}

sub _parse_next_res {
    my $self    = shift;
    my $s       = $self->settings;
    my $query = $self->_prepare_query();

    my $indexfile = $s->datapath . '/' . $s->database . '.idx';
    my $idx = Bio::Index::Fasta->new($indexfile);
    
    my %mapping = $self->_mapping;

    my $FH = $self->_output_fh;
    while (my $line = <$FH>) {
        chomp $line;
        my ( $pos, $sequence ) = split ":", $line;
        my $id = $mapping{$pos};

        my $seq_subject = $idx->fetch($id);

        my $seq_query = Bio::Seq->new(
            -display_id => "Query",
            -seq        => $query
        );

        # my $seq_subject = Bio::Seq->new(
        #    -display_id => "Subject",
        #    -seq        => $query
        # );

        my $alignment = undef;
        $alignment = $self->_get_alignment( $seq_query, $seq_subject )
            unless $s->no_alignments;

        my $res = Bio::Grep::SearchResult->new( $seq_subject,
            $s->upstream, $s->upstream + length($query),
            $alignment, $seq_subject->id, "" );
        # agrep does not support multiple queries yet    
        $res->query(Bio::Seq->new( -display_id => "Query", -seq => $s->query));
        return $res;
    }
    $self->_delete_output();
    return 0;
}

sub get_sequences {
    my ( $self, $seqid ) = @_;
    return $self->_get_sequences_from_bio_index($seqid);
}

sub available_sort_modes {
    return ();
}
1;    # Magic true value required at end of module
__END__


=head1 NAME

Bio::Grep::Backend::Agrep - Agrep back-end  


=head1 SYNOPSIS

  use Bio::Grep;
  
  my $sbe = Bio::Grep->new('Agrep');
  
  $sbe->settings->datapath('data');
  
  # generate a database. you have to do this only once. 
  $sbe->generate_database('ATH1.cdna', 'AGI Transcripts (- introns, + UTRs)');
  
  # search for the reverse complement and allow 2 mismatches 
  # Don't calculate Alignments with EMBOSS
  $sbe->search({
    query   => 'GAGCCCTT',
    reverse_complement => 1, 
    mismatches         => 2,
    no_alignments      => 1,
    database           => 'ATH1.cdna',
  });
  
  my @internal_ids;
  
  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     # print $res->alignment_string() . "\n\n";
     push @internal_ids, $res->sequence_id;
  }
  
  # get the complete sequences as Bio::SeqIO object
  my $seq_io = $sbe->get_sequences(\@internal_ids);

=head1 DESCRIPTION

B<Bio::Grep::Backend::Agrep> searches for a query with agrep. 

Note: L<Bio::Grep::Backend::Agrep> databases are compatible with
L<Bio::Grep::Backend::RE> databases.

=head1 METHODS

See L<Bio::Grep::Backend::BackendI> for other methods. 

=over 2

=item C<Bio::Grep::Backend::Agrep-E<gt>new()>

This function constructs an agrep back-end object

   my $sbe = Bio::Grep::Backend::Agrep->new();

=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

Available sortmodes in Agrep:

=over

   currently none.
   
=back

=back


=head1 DIAGNOSTICS

See L<Bio::Grep::Backend::BackendI> for other diagnostics. 

=over

=item C<Agrep call failed. Command was: ...> 

It was not possible to run agrep in function search(). Check the search
settings. Agrep returns also exit(1) whenever no hit is found! If you want to
reproduce the system() call, you can set the environment variable
C<BIOGREPDEBUG>. If this variable is set, then the temporary files won't get
deleted.  C<Bio::Root::SystemException>.

=back

=head1 SEE ALSO

L<Bio::Grep::Backend::BackendI>
L<Bio::Grep::SearchSettings>
L<Bio::SeqIO>
L<Bio::Index::Fasta>


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
