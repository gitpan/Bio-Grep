package Bio::Grep::Backends::Agrep;

use strict;
use warnings;

use Bio::Grep::Container::SearchSettings;
use Bio::Grep::Container::SearchResult;
use Bio::Grep::Backends::BackendI;

use base 'Bio::Grep::Backends::BackendI';

use File::Basename;
use IO::String;

use version; our $VERSION = qv('0.5.1');

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
    delete $all_features{QUERYFILE};
    delete $all_features{SHOWDESC};
    delete $all_features{QSPEEDUP};
    delete $all_features{REVCOM_DEFAULT};
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
    $mm = $s->mismatches if $s->mismatches > 0;

    # make insertions and deletions to expensive
    my $fuzzy = "-$mm -D" . ( $mm + 1 ) . " -I" . ( $mm + 1 );

    if ( $s->editdistance_isset ) {
        $fuzzy = "-" . $s->editdistance;
        $self->warn(
            "Ignoring mismatch settings because edit distance is set.")
            if $s->mismatches > 0;
    }

    $s->query_length( length($query) ) unless $s->query_length_isset;
    my $length = "-l " . $s->query_length;

    my $command =
        $self->_cat_path_filename( $s->execpath, 'agrep' ) . ' -i ' . $fuzzy
        . ' '
        . $query . ' '
        . $self->_cat_path_filename( $s->datapath, $s->database . '.dat' );

    if ( $ENV{BIOGREPDEBUG} ) {
        warn $command . "\n";
    }

    my $output = `$command`;
    my $cmd_ok = $self->_execute_command($command);

    # delete temporary files
    #unlink($tmp_query_file) if !$query_file;

    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  => "Agrep error: Query not valid.\n\nCommand was:\n\t$command"
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

sub generate_database_out_of_fastafile {
    my ( $self, $file, $description ) = @_;
    my ( $filename ) = fileparse($file);

    $self->_copy_fasta_file_and_create_nfo( $file, $filename, $description );

    $filename = $self->_cat_path_filename($self->settings->datapath,
        $filename);

    open my $DATFILE, '>',  "$filename.dat" 
        or $self->throw(
        -class => 'Bio::Root::FileOpenException',
        -text  => "Can't open $filename.dat for writing",
        -value => $!
        );
    
    open my $MAPFILE, '>', "$filename.map"
        or $self->throw(
        -class => 'Bio::Root::FileOpenException',
        -text  => "Can't open $filename.map for writing",
        -value => $!
        );

    my $in = Bio::SeqIO->new( -file => $filename, -format => 'Fasta' );
    my $i = 0;
    
    while ( my $seq = $in->next_seq() ) {
        print $MAPFILE $i . "\t" . $seq->id
            . "\n"
            or $self->throw(
            -class => 'Bio::Root::IOException',
            -text  => "Can't write to $filename.map",
            -value => $!
            );
        print $DATFILE $i . ":" . $seq->seq
            . "\n"
            or $self->throw(
            -class => 'Bio::Root::IOException',
            -text  => "Can't write to $filename.dat",
            -value => $!
            );
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

    open my $MAPFILE, '<', $mapfile 
        or $self->throw(
        -class => 'Bio::Root::FileOpenException',
        -text  => "Can't open $mapfile for reading",
        -value => $!
        );

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

        my $res = Bio::Grep::Container::SearchResult->new( $seq_subject,
            $s->upstream, $s->upstream + $s->query_length,
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

Bio::Grep::Backends::Agrep - Agrep back-end  


=head1 SYNOPSIS

  use Bio::Grep::Backends::Agrep;
  
  use Error qw(:try);
  
  # configure our search back-end, in this case Agrep
  my $sbe = Bio::Grep::Backends::Agrep->new();
  
  $sbe->settings->execpath('/usr/bin');
  $sbe->settings->tmppath('tmp');
  $sbe->settings->datapath('data_agrep');
  
  # generate a database. you have to do this only once. 
  $sbe->generate_database_out_of_fastafile('ATH1.cdna', 'AGI Transcripts (- introns, + UTRs)');
  
  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
  
  # take first available database in our test
  $sbe->settings->database($local_dbs[0]);
  
  my $seq = 'UGAACAGAAAGCUCAUGAGCC'; 
  
  # search for the reverse complement and allow 4 mismatches
  $sbe->settings->query($seq);
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(4);
  
  # things should be fast? then turn alignment calculation off
  $sbe->settings->no_alignments(1);
  
  # now try to search, agrep back-end will also throw an exception if no hit
  # was found
  try {
     $sbe->search();
  } catch Error::Simple with {
     my $E = shift;
     print  $E->{'-text'} . ' (' .  $E->{'-line'} . ")";
  };
  
  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";
  }


=head1 DESCRIPTION

B<Bio::Grep::Backends::Agrep> searches for a query with agrep. 


=head1 METHODS

See L<Bio::Grep::Backends::BackendI> for other methods. 

=over 2

=item C<Bio::Grep::Backends::Agrep-E<gt>new()>

This function constructs an agrep back-end object

   my $sbe = Bio::Grep::Backends::Agrep->new();

=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

Available sortmodes in Agrep:

=over

   currently none.
   
=back

=back


=head1 DIAGNOSTICS

See L<Bio::Grep::Backends::BackendI> for other diagnostics. 

=over

=item C<Bio::Root::SystemException>

It was not possible to run agrep in function C<search>. Check the search
settings. Agrep returns also exit(1) if no hit was found!

=item C<Bio::Root::FileOpenException>

=over 2

=item In function C<search>: It was not possible to open the database. Check
permissions and paths.

=item In function C<generate_database_out_of_fastafile>: It was not possible to
write the database to disk. Check permissions and paths.

=back


=item C<Bio::Root::IOException>

It was not possible to write the database in
C<generate_database_out_of_fastafile>. Check free diskspace.

=back


=head1 SEE ALSO

L<Bio::Grep::Backends::BackendI>
L<Bio::Grep::Container::SearchSettings>
L<Bio::SeqIO>
L<Bio::Index::Fasta>


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
