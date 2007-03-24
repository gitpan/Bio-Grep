package Bio::Grep::Backends::Hypa;

use strict;
use warnings;

use Bio::Grep::Container::SearchResult;
use Bio::Grep::Backends::BackendI;

use base 'Bio::Grep::Backends::BackendI';

use File::Temp qw/ tempfile tempdir /;
use File::Basename;

use version; our $VERSION = qv('0.5.0');

sub new {
    my $self = shift;
    $self = $self->SUPER::new;
    my %all_features = $self->features;
    delete $all_features{NATIVE_ALIGNMENTS};
    delete $all_features{EVALUE};
    delete $all_features{PERCENT_IDENTITY};
    delete $all_features{EDITDISTANCE};
    delete $all_features{SORT};
    delete $all_features{COMPLETE};
    delete $all_features{SHOWDESC};
    delete $all_features{QSPEEDUP};
    delete $all_features{HXDROP};
    delete $all_features{EXDROP};
    delete $all_features{REVCOM_DEFAULT};
    $self->features(%all_features);
    $self;
}

sub search {
    my ($self, $arg_ref) = @_;
    my $s    = $self->settings;
    $self->_check_search_settings($arg_ref);

    my $query   = $self->_prepare_query();

    my ( $tmp_fh, $tmp_query_file ) = tempfile(
        'spatter_XXXXXXXXXXXX',
        DIR    => $s->tmppath,
        SUFFIX => '.hypa'
    );

    $self->_make_hypa_query_file( $tmp_fh );

    my $online = '';
    $online = ' -online ' if ( $s->online_isset && $s->online );

    my $upstream = '';
    $upstream = ' -showleft ' . $s->upstream . ' ' if $s->upstream > 0;
    my $downstream = '';
    $downstream = ' -showright ' . $s->downstream . ' ' if $s->downstream > 0;
    #set query_length automatically
    my $auto_query_length = 0;
    if (!defined $s->query_length) {
        $s->query_length( length($s->query) );
        $auto_query_length = 1;
    }    

    my $command =
        $self->_cat_path_filename( $s->execpath, 'hypa' )
        . " -showformat -nohitcount -delimiter ';' -q "
        . $tmp_query_file
        . $online
        . $upstream
        . $downstream
        . ' -index '
        . $self->_cat_path_filename( $s->datapath, $s->database );

    #warn $command . " \n ";
    my $cmd_ok = $self->_execute_command($command);

    # delete temporary files
    unlink($tmp_query_file);

    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  => "Hypa: Maybe query not valid. Command was:\n\t$command"
        )
        if !$cmd_ok;

    #$self->_output_parser($output);
    #$self->_filter_results();
    $self->_skip_header();
    $self->_prepare_results;
    $self->settings->query_length_reset if $auto_query_length;
    return 1;
}

sub _make_hypa_query_file {
    my ( $self, $file_fh ) = @_;

    my $editops = "";
    my $cost    = "";
    my $query   = $self->_prepare_query();
    my $gu_mm   = 1;
    $gu_mm = $self->settings->gumismatches
        if $self->settings->gumismatches_isset;
    my $mismatches = 0;
    $mismatches = $self->settings->mismatches
        if $self->settings->mismatches_isset;
    my $deletions = 0;
    $deletions = $self->settings->deletions
        if $self->settings->deletions_isset;
    my $insertions = 0;
    $insertions = $self->settings->insertions
        if $self->settings->insertions_isset;

    if ( $gu_mm != 0.5 && $gu_mm != 1 && $gu_mm != 0 ) {
        $self->warn("Valid values for GU mismatches are: 0, 0.5 and 1");
    }

    my $gu1 = 'G->T';
    my $gu2 = 'T->G';

    $self->settings->reverse_complement(0)
        unless $self->settings->reverse_complement_isset;
    if ( $self->settings->reverse_complement ) {
        $gu1 = 'A->G';
        $gu2 = 'C->T';
    }

    # what shall we do with GU?
    if ( $gu_mm == 1 ) {    # normal mode,  a mismatch
        $editops = "[$mismatches, $deletions, $insertions]";
    }
    elsif ( $gu_mm == 0 ) {    # not a mismatch
        $cost = <<EOT;
// match 0 * 2
cost(T->T) := 0;
cost(A->A) := 0;
cost(C->C) := 0;
cost(G->G) := 0;
cost($gu1) := 0;
cost($gu2) := 0;

// insertions, deletions, mismatches count 1 * 2
cost(.->) := 2;
cost(->.) := 2;
cost(.->T) := 2;
cost(.->C) := 2;
cost(.->G) := 2;
cost(.->A) := 2;

EOT
        $editops = "[cost, "
            . $mismatches * 2 . ","
            . $deletions * 2 . ","
            . $insertions * 2 . "]";
    }
    else {
        $cost = <<EOT;
// match 0 * 2
cost(T->T) := 0;
cost(A->A) := 0;
cost(C->C) := 0;
cost(G->G) := 0;

// gu 1/2 * 2 
cost($gu1) := 1;
cost($gu2) := 1;

// insertions, deletions, mismatches count 1 * 2
cost(.->) := 2;
cost(->.) := 2;
cost(.->T) := 2;
cost(.->C) := 2;
cost(.->G) := 2;
cost(.->A) := 2;

EOT
        $editops = "[cost, "
            . $mismatches * 2 . ","
            . $deletions * 2 . ","
            . $insertions * 2 . "]";
    }

    $query = "($query)";

    my $hypaquery = <<EOT;
   $cost
   query = $query$editops;

EOT

    print $file_fh $hypaquery;
    close $file_fh;
}

sub get_databases {
    my $self = shift;
    return $self->_get_databases('.rev.suf');
}

sub generate_database_out_of_fastafile {
    my $self        = shift;
    my $file        = shift;
    my $description = shift;
    my ( $filename, $oldpath ) = fileparse($file);

    $self->_copy_fasta_file_and_create_nfo( $file, $filename, $description );

    my $filepath = $self->_cat_path_filename($self->settings->datapath,
        $filename);
    my $alphabet = $self->_guess_alphabet_of_file($file);
    my $alphabet_specific_arguments = '';
    #warn $alphabet; 
    if ($alphabet eq 'protein') {
        $alphabet_specific_arguments = ' -protein ';
    }
    elsif ($alphabet eq 'dna') {
        $alphabet_specific_arguments = ' -dna ';
    }
    else {
        $self->throw(-class => 'Bio::Root::BadParameter',
                     -text  => 'unsupported alphabet of file',
                     -value => $alphabet,);   
    }


    my $command =
        $self->_cat_path_filename( $self->settings->execpath, 'mkaffix.sh' )
        . ' '
        . $filename . ' '
        . $filename
        . $alphabet_specific_arguments;    
                      #print $command;
    my $output_dir = $self->settings->datapath;
    system(qq{ cd $output_dir ; exec $command } );
    
    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  =>
            "Hypa error: Cannot generate suffix array. Command was:\n\t$command"
        )
        if ($?);

    # create database from directory of fasta files
    my $idx = Bio::Index::Fasta->new(
        -filename   => $filepath . '.idx',
        -write_flag => 1
    );
    $idx->make_index( ($filepath) );
    return $filename;    
}

sub _skip_header {
    my ( $self ) = @_;
    my $FH = $self->_output_fh;
    while (my $line = <$FH>) {
        chomp $line;
        
        # skip everything before first line
        if ($line =~ /^\s*$/) {
            return;
        }
    }
}

sub _parse_next_res {
    my $self    = shift;
    my $s       = $self->settings;
    my $query   = $self->_prepare_query();

    # temp variables. for parsing only
    my $subject = '';
    my $tmp_seq;
    my $upstream   = '';
    my $match      = '';
    my $downstream = '';

    my $FH = $self->_output_fh;
    while (<$FH>) {
        chomp;
        last if (/^Matches for pattern/);
        if (/Time/) { print "$_\n"; last; }

        # skip everything before first line
        if (/^\s*$/) {
            next;
        }
        if (/^>(.*)\:\[(\d+)\,(\d+)\]/) {
            $tmp_seq = Bio::LocatableSeq->new(
                -id  => $1,
                -seq => 'AAA'
                , # this is just a dummy, we will overwrite this when we store the hit
                -start => $2,
                -end   => $3
            );
            $subject = '';
            next;
        }
        my @parts = split (';');

        $upstream   = '';
        $match      = '';
        $downstream = '';
        foreach my $part (@parts) {
            $upstream   = $part if $part =~ /^\^/;
            $match      = $part if $part =~ /^[^\^\$]/;
            $downstream = $part if $part =~ /^\$/;
        }
        $upstream   =~ s/^\^//;
        $downstream =~ s/^\$//;
        next
            if $match eq
            '';   # this should not be, but simply ignore this case, otherwise
                  # parsing would fail
                  # first newline after hit, so now store hit.
        my $seq_query = Bio::Seq->new(
            -display_id => 'Query',
            -seq        => $query
        );

        # align only match
        my $seq_subject = Bio::Seq->new(
            -display_id => $tmp_seq->id,
            -seq        => $match
        );
        $tmp_seq->seq( $upstream . $match . $downstream );
        my $alignment = undef;
        $alignment = $self->_get_alignment( $seq_query, $seq_subject )
            unless $s->no_alignments;
        my $res = $self->_filter_result(
            Bio::Grep::Container::SearchResult->new( $tmp_seq,
            length($upstream), length($upstream) + length($match),
            $alignment, $tmp_seq->id, '' )
            );
        if ($res) {    
            $res->query(Bio::Seq->new( -display_id => "Query", -seq => $s->query));
            return $res;
        }    
    }
    $self->_delete_output();
    return 0;
}

sub get_sequences {
    my ( $self, $seqid ) = @_;
    return $self->_get_sequences_from_bio_index($seqid);
}

sub available_sort_modes {
    my ($self) = @_;

    # get sort modes from superclass
    return ( $self->SUPER::available_sort_modes() );
}
1;    # Magic true value required at end of module
__END__


=head1 NAME

Bio::Grep::Backends::Hypa - HyPa back-end


=head1 SYNOPSIS

  use Bio::Grep::Backends::Hypa;
 
  use Bio::Root::Exception;
  use Error qw(:try);
 
  # construct the Hypa back-end	
  my $sbe = Bio::Grep::Backends::Hypa->new();
 
  $sbe->settings->tmppath('tmp');
  $sbe->settings->datapath('data');
  
  # generate a Hypa suffix array. you have to do this only once.
  $sbe->generate_database_out_of_fastafile('ATH1.cdna', 'AGI Transcripts (- introns, + UTRs)');
 
  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
 
  # take first available database in our test
  $sbe->settings->database($local_dbs[0]);
 
  my $seq = 'UGAACAGAAAGCUCAUGAGCC'; 
 
  # search for the reverse complement and allow 4 mismatches
  # display 5 bases upstream and downstream of the match
  $sbe->settings->query($seq);
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(4);
  $sbe->settings->upstream(5);
  $sbe->settings->downstream(5);
  
  # hypa does not generate alignments, but EMBOSS will automatically calculate some for
  # you. If things should be fast, turn this off	
  # $sbe->settings->no_alignments(1);
 
  # with the features hash, you can check if the back-end supports
  # a special feature. Allow wobble pairs (GU counts 0 or 0.5)	
  if (defined($sbe->features->{GUMISMATCHES})) {
 	  $sbe->settings->gumismatches(0.5); 		
  } else {
     print "\nBack-end does not support wobble pairs\n";
  }	
 
  # now try to search. we use bioperls Exceptions
  try {
     $sbe->search();
  } catch Bio::Root::SystemException with {
     my $E = shift;
     print STDERR 'Back-end call failed: ' . 	
     $E->{'-text'} . ' (' .  $E->{'-line'} . ")\n";
     exit(1);	
  } otherwise {        
 	my $E = shift;
     	print STDERR "An unexpected exception occurred: \n$E";
     exit(1);	
 
  };
 
  # output all informations we have!
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";
  }
   
=head1 DESCRIPTION

B<Bio::Grep::Backends::Hypa> searches for a query in a Hypa suffix array. 

NOTE 1: Hypa can not calculate alignments. But because we have the exact position
of the match, the alignment calculation shouldn't be too slow (In agrep, we align query 
and sequence, here query and approximate match).

NOTE 2: -online is available. But it is not recommended.

=head1 METHODS

See L<Bio::Grep::Backends::BackendI> for other methods. 

=over 2

=item Bio::Grep::Backends::Hypa-E<gt>new()

This function constructs a Hypa back-end object

   my $sbe = Bio::Grep::Backends::Hypa->new();

=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

   $sbe->sort('ga');

Available sortmodes in Hypa:

=over

            ga  : 'ascending order of dG'
            gd  : 'descending order of dG'

=back

Note that 'ga' and 'gd' require that search results have dG set. 
L<Bio::Grep::RNA> ships with filters for free energy calculation. Also note that
these two sort options require that we load all results in memory.

=back


=head1 BUGS AND LIMITATIONS

Hypa currently does not return the correct downstream and upstream regions if
the requested region is larger then available (if hit is near 5' or 3' end).

Hypa does not generate alignments so we use EMBOSS for that. This is slow and
produces a lot of IO.

Please report any bugs or feature requests to
C<bug-bio-grep@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>. 


=head1 SEE ALSO

L<Bio::Grep::Backends::BackendI>
L<Bio::Grep::Container::SearchSettings>
L<Bio::SeqIO>


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
