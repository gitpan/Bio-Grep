package Bio::Grep::Backend::Vmatch;

use strict;
use warnings;

use Fatal qw(open close);

use Bio::Grep::SearchResult;
use Bio::Grep::Backend::BackendI;

use base 'Bio::Grep::Backend::BackendI';

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use IO::String;

use version; our $VERSION = qv('0.9.1');

sub new {
    my $self = shift;
    $self = $self->SUPER::new;
    my %all_features = $self->features;
    delete $all_features{GUMISMATCHES};
    delete $all_features{DELETIONS};
    delete $all_features{INSERTIONS};
    delete $all_features{REVCOM_DEFAULT};
    $self->features(%all_features);
    $self;
}

sub search {
    my ( $self, $arg_ref ) = @_;
    my $s = $self->settings;
    $self->_check_search_settings($arg_ref);

    my ( $query, $query_file, $tmp_query_file )
        = $self->_create_tmp_query_file();

    if ( ( $s->upstream > 0 || $s->downstream > 0 ) && $s->showdesc_isset ) {
        $self->throw(
            -class => 'Bio::Root::BadParameter',
            -text  => "You can't use showdesc() with upstream or downstream.",
        );
    }
    if ( $query_file && !$s->complete && !$s->query_length_isset ) {
        $self->throw(
            -class => 'Bio::Root::BadParameter',
            -text  => 'You have to specify complete or querylength. See '
                . 'the flags -complete and -l in the Vmatch documentation.',
        );
    }

    # now generate the command string
    my $fuzzy = '';
    if ( $s->mismatches_isset && $s->mismatches > 0 ) {
        $fuzzy = ' -h ' . $s->mismatches . ' ';
    }
    if ( $s->editdistance_isset && $s->editdistance > 0 ) {
        $fuzzy = ' -e ' . $s->editdistance . ' ';
    }
    my $online = '';
    $online = ' -online ' if ( $s->online_isset && $s->online );
    my $auto_query_length = 0;
    if ( !defined $s->query_length && !$s->complete && !$query_file ) {
        $s->query_length( length($query) );
        $auto_query_length = 1;
    }
    my $length = '';
    $length = ' -l ' . $s->query_length if defined $s->query_length;
    my $complete = '';
    $complete = ' -complete ' if ( $s->complete_isset && $s->complete );

    my $sort = '';

    # dg sorting done by BackendI, not vmatch
    if ( $s->sort_isset && substr( $s->sort, 0, 1 ) ne 'g' ) {
        $sort = ' -sort ' . $s->sort;
    }
    my $maxhits = '';
    $maxhits = ' -best ' . $s->maxhits . ' ' if $s->maxhits_isset;

    my $showdesc = '';
    $showdesc = ' -showdesc ' . $s->showdesc . ' ' if $s->showdesc_isset;

    my $qspeedup = '';
    if ( $s->qspeedup_isset ) {
        $self->throw(
            -class => 'Bio::Root::BadParameter',
            -text  => "You can't combine qspeedup and complete.",
        ) if $s->complete_isset;

        $qspeedup = ' -qspeedup ' . $s->qspeedup . ' ';

    }

    my $hxdrop = '';
    $hxdrop = ' -hxdrop ' . $s->hxdrop . ' ' if $s->hxdrop_isset;

    my $exdrop = '';
    $exdrop = ' -exdrop ' . $s->exdrop . ' ' if $s->exdrop_isset;

    my $pflag = '';
    $pflag = ' -p ' if ( $s->query_file && $s->reverse_complement );

    if ( $s->direct_and_rev_com ) {
        $pflag = ' -p -d ';
    }
    my $command = $self->_cat_path_filename( $s->execpath, 'vmatch' ) . ' -q '
        . $tmp_query_file
        . $complete
        . $sort
        . $fuzzy
        . $maxhits
        . $qspeedup
        . $showdesc
        . $length . ' -s '
        . $online
        . $hxdrop
        . $exdrop
        . $pflag
        . $self->_cat_path_filename( $s->datapath, $s->database );

    if ( $ENV{BIOGREPDEBUG} ) {
        warn $command . "\n";
    }

    my $cmd_ok = $self->_execute_command($command);

    # delete temporary files
    #unlink($tmp_query_file) if !$query_file;

    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  => "Vmatch call failed. Command was:\n\t$command"
    ) if !$cmd_ok;

    if ( $s->showdesc_isset ) {
        my %query_desc_lookup;
        foreach my $query_seq ( @{ $self->_query_seqs } ) {

            # simulate how this sequence would look in vmatch output
            my $query_desc = $query_seq->id;
            $query_desc .= ' ' . $query_seq->desc
                if ( $query_seq->desc
                && length($query_desc) > 0 );
            $query_desc = substr $query_desc, 0, $s->showdesc;
            $query_desc =~ s{ }{_}g;
            $query_desc_lookup{$query_desc} = $query_seq;
        }
        $self->_mapping(%query_desc_lookup);
    }
    $self->settings->query_length_reset if $auto_query_length;
    $self->_prepare_results;
    return 1;
}

sub get_databases {
    my $self = shift;
    return $self->_get_databases('.al1');
}

sub generate_database {
    my ( $self, @args ) = @_;
    my %args = $self->_prepare_generate_database(@args);

    if (defined $args{skip}) {
        return 0;
    }   

    my $alphabet = $self->_guess_alphabet_of_file( $args{filename} );
    my $alphabet_specific_arguments = '';

    #warn $alphabet;
    #
    my $verbose = '';
    $verbose = ' -v ' if defined $args{verbose};

    my $pl = ' -pl ';
    $pl .= $args{prefix_length} . ' ' if defined $args{prefix_length};

    if ( $alphabet eq 'protein' ) {
        $alphabet_specific_arguments = ' -protein ';
    }
    elsif ( $alphabet eq 'dna' ) {
        $alphabet_specific_arguments = ' -dna ';
    }
    else {
        $self->throw(
            -class => 'Bio::Root::BadParameter',
            -text  => 'Unsupported alphabet of file.',
            -value => $alphabet,
        );
    }
    my $command
        = $self->_cat_path_filename( $self->settings->execpath, 'mkvtree' )
        . ' -db '
        . $args{basefilename}
        . $alphabet_specific_arguments
        . $pl
        . ' -allout '
        . $verbose;

    if ( $ENV{BIOGREPDEBUG} ) {
        warn $command . "\n";
    }
    my $output_dir = $self->settings->datapath;
    system(qq{ cd $output_dir ; exec $command });
    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text =>
            "mkvtree call failed. Cannot generate suffix array. Command was:\n\t$command"
    ) if ($?);
    return 1;
}

sub _parse_next_res {
    my $self                = shift;
    my @query_seqs          = $self->_query_seqs;
    my $s                   = $self->settings;
    my @results             = ();
    my $alignment_in_output = -1;
    my $skip_next_alignment = 0;
    my $skipped_lines       = 0;
    my $subject;
    my $tmp_aln;

    my ( $command, $output );

    my $FH = $self->_output_fh;
LINE:
    while ( my $line = <$FH> ) {
        chomp $line;
        $line =~ s/\s+/ /g;
        if ( $line !~ /\s/ ) {
            $skipped_lines++;
            if ( $skipped_lines == 2 ) {
                $results[-1]->alignment($tmp_aln);
                my $real_subject = $subject;

                # remove gaps out of alignment
                $real_subject =~ s{-}{}g;
                $results[-1]->sequence->seq($real_subject)
                    if $s->showdesc_isset;

                #            $results[-1]->query->seq($query);
                my $res = $self->_filter_result( $results[-1] );
                return $res if $res;
                $alignment_in_output = -1;
                next LINE;
            }
            else {
                next LINE;
            }
        }
        else {
            $skipped_lines = 0;
        }
        my @fields = split ' ', $line;
        $alignment_in_output = 0
            if ( $line =~ /^Sbjct:/ && !$skip_next_alignment );

        next unless ( $fields[0] =~ /^\d+$/ || $alignment_in_output >= 0 );
        $skip_next_alignment = 0;

        if ( $line =~ /^Sbjct: (.*)$/ ) {
            $subject = $1;
        }
        if ( $line =~ /^Query: (.*)$/ ) {
            my $query = $1;
            $query =~ s/\s+(\d+)\s*$//;
            my $query_pos = $1;
            $subject =~ s/\s+(\d+)\s*$//;
            my $subject_pos = $1;
            if ( !$tmp_aln->no_sequences ) {
                $tmp_aln->add_seq(
                    Bio::LocatableSeq->new(
                        -id    => 'Subject',
                        -seq   => $subject,
                        -start => ( $subject_pos - length($subject) ) + 1,
                        -end   => $subject_pos
                    )
                );
                $tmp_aln->add_seq(
                    Bio::LocatableSeq->new(
                        -id    => "Query",
                        -seq   => $query,
                        -start => ( $query_pos - length($query) ) + 1,
                        -end   => $query_pos
                    )
                );
            }
            else {
                my $s1 = $tmp_aln->get_seq_by_pos(1);
                my $s2 = $tmp_aln->get_seq_by_pos(2);
                $s1->end($subject_pos);
                $s1->seq( $s1->seq . $subject );
                $s2->end($query_pos);
                $s2->seq( $s2->seq . $query );
                $tmp_aln = new Bio::SimpleAlign( -source => "VMATCH" );
                $tmp_aln->add_seq($s1);
                $tmp_aln->add_seq($s2);
            }
            next LINE;
        }
        next LINE if $alignment_in_output >= 0;

        $tmp_aln = new Bio::SimpleAlign( -source => "VMATCH" );

        # make taint mode happy
        ( $fields[0] ) = $fields[0] =~ /(\d+)/;
        ( $fields[2] ) = $fields[2] =~ /(\d+)/;

        if ( $s->showdesc_isset ) {

            # not numerical with showdesc on
            # don't worry about this unclean untainting, we don't use that
            # description in dangerous ways
            ( $fields[1] ) = $fields[1] =~ /(.+)/;
            ( $fields[5] ) = $fields[5] =~ /(.+)/;
        }
        else {
            ( $fields[1] ) = $fields[1] =~ /(\d+)/;
            ( $fields[5] ) = $fields[5] =~ /(\d+)/;
        }

        #     warn Data::Dumper->Dump([ %descs ]);

        my $fasta;
        my $upstream        = $s->upstream;
        my $internal_seq_id = $fields[1];

        if ( !$s->showdesc_isset ) {
            my $start = $fields[2] - $s->upstream;

            # maybe the defined upstream region is larger than available
            # so check this and store in local variables
            if ( $start < 0 ) {
                $upstream = $upstream + $start;
                $start    = 0;
            }
            my $length = $upstream + $fields[0] + $s->downstream;
            $command
                = $self->_cat_path_filename( $s->execpath, 'vsubseqselect' )
                . " -seq $length $fields[1] "
                . $start . ' '
                . $self->_cat_path_filename( $s->datapath, $s->database );
            $output = `$command`;
            if ( $ENV{BIOGREPDEBUG} ) {
                warn $command . "\n";
            }
            my $stringio = IO::String->new($output);
            my $in       = Bio::SeqIO->new(
                '-fh'     => $stringio,
                '-format' => 'fasta'
            );
            $fasta = $in->next_seq();
        }
        else {
            my ( $seq_id, $seq_desc )
                = $fields[1] =~ m{\A (.+?) _ (.*) \z}xms;
            $seq_id   = $fields[1] if !defined $seq_id;
            $seq_desc = ''         if !defined $seq_desc;
            $seq_desc =~ s{_}{ }g;
            $fasta = Bio::Seq->new( -id => $seq_id, -desc => $seq_desc );
            $internal_seq_id = $seq_id;
        }

        my $query;
        if ( $s->showdesc_isset ) {
            $query = $self->_mapping->{ $fields[5] };
        }
        else {
            $query = $query_seqs[ $fields[5] ];
        }

        my $result = Bio::Grep::SearchResult->new(
            {   sequence         => $fasta,
                begin            => $upstream,
                end              => $upstream + $fields[0],
                alignment        => Bio::SimpleAlign->new(),
                sequence_id      => $internal_seq_id,
                remark           => '',
                evalue           => $fields[8],
                percent_identity => $fields[10],
                query            => $query,
            }
        );
        push( @results, $result );
    }
    $self->_delete_output();
    return 0;
}

sub get_sequences {
    my ( $self, $seqid ) = @_;
    my $s = $self->settings;
    $self->is_arrayref_of_size( $seqid, 1 );
    $self->_check_search_settings();
    my ( $tmp_fh, $tmpfile );

    my $seq_query = '';

    if ( @{$seqid}[0] =~ m{\A \d+ \z}xms ) {
        ( $tmp_fh, $tmpfile )
            = tempfile( 'vseqselect_XXXXXXXXXXXXX', DIR => $s->tmppath );

        foreach ( @{$seqid} ) {
            print $tmp_fh $_ . " \n ";
        }
        close $tmp_fh;
        $seq_query = ' -seqnum ' . $tmpfile;
    }
    else {
        my $seq_desc = $self->is_sentence( @{$seqid}[0] );
        $seq_query = ' -matchdesc "' . $seq_desc . '"';
    }

    my $command = $self->_cat_path_filename( $s->execpath, 'vseqselect' )
        . $seq_query . ' '
        . $self->_cat_path_filename( $s->datapath, $s->database );

    if ( $ENV{BIOGREPDEBUG} ) {
        warn $command . "\n";
    }

    my $output = `$command`;
    if ( $? && $output !~ m{\A \> }xms ) {
        $self->throw(
            -class => 'Bio::Root::SystemException',
            -text =>
                "vseqselect call failed. Cannot fetch sequences. Command was:\n\t$command\n$output"
        );
    }
    unlink($tmpfile) if $tmpfile;
    my $stringio = IO::String->new($output);
    my $out      = Bio::SeqIO->new(
        '-fh'     => $stringio,
        '-format' => 'fasta'
    );
    return $out;
}

sub available_sort_modes {
    my ($self) = @_;
    return (
        $self->SUPER::available_sort_modes(),
        la  => 'ascending order of length',
        ld  => 'descending order of length',
        ia  => 'ascending order of first position',
        id  => 'descending order of first position',
        ja  => 'ascending order of second position',
        jd  => 'descending order of second position',
        ea  => 'ascending order of Evalue',
        ed  => 'descending order of Evalue',
        sa  => 'ascending order of score',
        sd  => 'descending order of score',
        ida => 'ascending order of identity',
        idd => 'descending order of identity'
    );
}
1;    # Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep::Backend::Vmatch - Vmatch back-end  

=head1 SYNOPSIS

  use Bio::Grep;
  
  my $sbe = Bio::Grep->new('Vmatch');
  
  # generate a Vmatch suffix array. you have to do this only once.
  $sbe->generate_database({ 
    file          => 'ATH1.cdna', 
    description   => 'AGI Transcripts',
    datapath      => 'data',
    prefix_length => 3,
  });
 
  # search for the reverse complement and allow 4 mismatches
  # parse the description (max. 100 chars) directly out of the
  # Vmatch output instead of calling vsubseqselect for every
  # search result

  $sbe->search({
    query   => 'UGAACAGAAAGCUCAUGAGCC',
    reverse_complement => 1,
    mismatches         => 4,
    showdesc           => 100,
    database           => 'ATH1.cdna',
  });

  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";

     # sequence_id now contains the gene id (e.g. At1g1234),
     # not the Vmatch internal id 
     # To retrieve the complete sequences, one has to
     # call get_sequences for every gene id
     my $seq_io = $sbe->get_sequences([$res->sequence_id]);
     my $sequence = $seq_io->next_seq;
  }
  
  # for retrieving up- and downstream regions,
  # Vmatch internal sequence ids are required
  # (no showdesc possible)

  $sbe->search({
    query   => 'AGAGCCCT',
    reverse_complement => 1,
    mismatches         => 1,
    upstream           => 30,
    downstream         => 30,
  });
 
  my @internal_ids;
  while ( my $res = $sbe->next_res ) {
    # vsubseqselect is called now for every result ...
    push @internal_ids, $res->sequence_id;
  }

  # ... but one can retrieve all complete sequences with
  # just one call of vseqselect
  my $seq_io = $sbe->get_sequences(\@internal_ids);
  
=head1 DESCRIPTION

B<Bio::Grep::Backend::Vmatch> searches for a query in a C<Vmatch> suffix array. 

=head1 METHODS

See L<Bio::Grep::Backend::BackendI> for inherited methods. 

=over 2

=item C<Bio::Grep::Backend::Vmatch-E<gt>new()>

This method constructs a C<Vmatch> back-end object and should not used 
directly. Rather, a back-end should be constructed by the main class 
L<Bio::Grep>:

  my $sbe = Bio::Grep->new('Vmatch');

=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

   $sbe->sort('ga');

Available sortmodes in C<Vmatch>:

=over

            ga  : 'ascending order of dG'
            gd  : 'descending order of dG'
            la  : 'ascending order of length'
            ld  : 'descending order of length'
            ia  : 'ascending order of first position'
            id  : 'descending order of first position'
            ja  : 'ascending order of second position'
            jd  : 'descending order of second position'
            ea  : 'ascending order of Evalue'
            ed  : 'descending order of Evalue'
            sa  : 'ascending order of score'
            sd  : 'descending order of score'
            ida : 'ascending order of identity'
            idd : 'descending order of identity'


=back

Note that 'ga' and 'gd' require that search results have dG set. 
L<Bio::Grep::RNA> ships with filters for free energy calculation.  Also note 
that these two sort options require that we load all results in memory.

=item C<$sbe-E<gt>get_sequences()>

Takes as argument an array reference. If first array element is an integer, 
then this method assumes that the specified sequence ids are C<Vmatch> internal 
ids. Otherwise it will take the first array element as query.

    # get sequences 0,2 and 4 out of suffix array
    $sbe->get_sequences([0,2,4]);

    # get sequences that start with At1g1
    $sbe->get_sequences(['At1g1', 'ignored']);

The internal ids are stored in C<$res-E<gt>sequence_id>. If you have specified
C<showdesc>, then C<sequence_id> will contain the gene id (e.g. At1g1234),
NOT the C<Vmatch> internal id.

=back

=head1 DIAGNOSTICS

See L<Bio::Grep::Backend::BackendI> for other diagnostics. 

=over

=item C<mkvtree call failed. Cannot generate suffix array. Command was: ...>. 

It was not possible to generate a suffix array in generate_database().
Check permissions and paths. C<Bio::Root::SystemException>.

=item C<Unsupported alphabet of file.>

The method generate_database() could not determine the
alphabet (DNA or Protein) of the specified Fasta file. C<Bio::Root::BadParameter>


=item C<Vmatch call failed. Command was: ... > 

It was not possible to run C<Vmatch> in function search(). Check the search
settings. If the number of mismatches is to high, try C<online>. 
C<Bio::Root::SystemException>.

=item C<vseqselect call failed. Cannot fetch sequences. Command was: ...> 

It was not possible to get some sequences out of the suffix array in 
get_sequences(). Check sequence ids. C<Bio::Root::SystemException>.

=item C<You can't combine qspeedup and complete.>

The C<Vmatch> parameters C<-complete> and C<-qspeedup> cannot combined. See 
the C<Vmatch> documentation. C<Bio::Root::BadParameter>.

=item C<You can't use showdesc() with upstream or downstream.>

We need the tool C<vsubseqselect> of the C<Vmatch> package for the upstream and 
downstream regions. This tool requires as parameter an internal C<Vmatch>
sequence id, which is not shown in the C<Vmatch> output when C<showdesc> is on. 
C<Bio::Root::BadParameter>.

=item C<You have to specify complete or querylength. ...'>

The C<Vmatch> parameters C<-complete> and C<-l> cannot combined. See the
C<Vmatch> documentation. C<Bio::Root::BadParameter>.

=back

=head1 SEE ALSO

L<Bio::Grep::Backend::BackendI>
L<Bio::Grep::SearchSettings>
L<Bio::SeqIO>

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
