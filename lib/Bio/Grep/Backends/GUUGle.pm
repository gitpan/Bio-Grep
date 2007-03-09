package Bio::Grep::Backends::GUUGle;

use strict;
use warnings;

use Bio::Grep::Container::SearchResult;
use Bio::Grep::Backends::BackendI;

use base 'Bio::Grep::Backends::BackendI';

use File::Basename;

use Data::Dumper;

use version; our $VERSION = qv('0.2.0');

sub new {
    my $self = shift;
    $self = $self->SUPER::new;
    my %all_features = $self->features;
    delete $all_features{EVALUE};
    delete $all_features{PERCENT_IDENTITY};
    delete $all_features{DELETIONS};
    delete $all_features{INSERTIONS};
    delete $all_features{EDITDISTANCE};
    delete $all_features{ONLINE};
    delete $all_features{SORT};
    delete $all_features{COMPLETE};
    delete $all_features{SHOWDESC};
    delete $all_features{QSPEEDUP};
    delete $all_features{PROTEINS};
    $self->settings->gumismatches(0);
    $self->features(%all_features);
    $self;
}

sub search {
    my ($self, $arg_ref) = @_;
    my $s    = $self->settings;
    $self->_check_search_settings($arg_ref);

    my ( $query, $query_file, $tmp_query_file) =
    $self->_create_tmp_query_file();

    # now generate the command string
    my $eflag = '';
    if ($s->upstream > 0 || $s->downstream > 0) {
        $self->throw( 
             -class => 'Bio::Root::BadParameter',
             -text => 'Upstream and downstream must be the same in GUUGle.',
             -value => 'Upstream: ' . $s->upstream . ' Downstream: ' .
             $s->downstream,
             ) if $s->upstream != $s->downstream;
        $eflag = ' -e ' . $s->upstream . ' ';     
    }
    
    if ($s->gumismatches > 0) {
        $self->warn('GUUGle counts GU always as no mismatch. ' .
            'Set gumismatches(0) to hide this warning.'
        );
    }   
    if (!$s->reverse_complement && $s->query_file) {
        $self->throw( 
             -class => 'Bio::Root::BadParameter',
             -text => 'GUUGle searches only for the reverse '. 
             'complement. ',
             -value => $s->reverse_complement,
         );     
    }   

    my $lflag = '';
    $lflag = ' -l ' . $s->maxhits . ' ' if $s->maxhits_isset;

    #set query_length automatically
    my $auto_query_length = 0;
    if (!defined $s->query_length) {
        if ($s->query) {
            $s->query_length( length($s->query) );
            $auto_query_length = 1;
        }
        else {
            $self->throw( 
                -class => 'Bio::Root::BadParameter',
                -text => 'settings->query_length not set. See -d flag in the '. 
                'GUUGle documentation. ',
         );     
        }    
    }    

    my $command =
        $self->_cat_path_filename( $s->execpath, 'guugle' )
        . $eflag
        . $lflag
        . ' -d ' . $s->query_length . ' '
        . $self->_cat_path_filename( $s->datapath, $s->database ) . ' '
        . $tmp_query_file
        ;

    if ( $ENV{BIOGREPDEBUG} ) {
        warn $command . "\n";
    }

    my $cmd_ok = $self->_execute_command($command);

    # delete temporary files
    #unlink($tmp_query_file) if !$query_file;

    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  => "GUUGle: Maybe query not valid. Command was:\n\t$command"
        )
        if !$cmd_ok;

    $self->_skip_header();
    $self->_prepare_results;
    my %query_desc_lookup;
    foreach my $query_seq (@{ $self->_query_seqs }) {
        # simulate how this sequence would look in guugle output
        my $query_desc = $query_seq->id;
        $query_desc .= ' ' . $query_seq->desc if length $query_seq->desc;
        $query_desc_lookup{$query_desc} = $query_seq;
    }
    $self->_mapping(%query_desc_lookup);
    $self->settings->query_length_reset if $auto_query_length;
    return 1;
}

sub get_databases {
    my $self = shift;
    return $self->_get_databases('.al1');
}

sub generate_database_out_of_fastafile {
    my ( $self, $file, $description ) = @_;
    my ( $filename ) = fileparse($file);

    $self->_copy_fasta_file_and_create_nfo( $file, $filename, $description );

    my $newfilename = $self->_cat_path_filename($self->settings->datapath,
        $filename);

    $self->_create_index_and_alphabet_file( $newfilename );
    return $filename;     
}

sub _skip_header {
    my ( $self ) = @_;
    my $FH = $self->_output_fh;
    my $blank_lines = 0;
    while (my $line = <$FH>) {
        chomp $line;
        # skip everything before first line
        if ($line =~ /^\s*$/) {
            $blank_lines++;
            return if $blank_lines == 2;
        }
    }
}

sub _rnas_match {
    my ( $self, $s1, $s2) = @_;
    return 1 if $s1 eq $s2;
    return 0 if length($s1) != length($s2);
    # now check for wobble pairs
    for my $i (0 .. (length($s1) -1)) {
        my $a = substr $s1, $i, 1;
        my $b = substr $s2, $i, 1;
        if ($a ne $b && join('',$a,$b) ne 'uc' && join('',$a,$b) ne 'ga') {
            return 0;
        }    
    }    
    return 1;
}    

sub _parse_next_res {
    my $self    = shift;
    my $s       = $self->settings;
    my @query_seqs = $self->_query_seqs;

    for my $i ( 0 .. $#query_seqs ) {

    }    
    # temp variables. for parsing only
    my ($subject_id,$subject_desc, $subject_pos, $subject_seq, 
        $query_desc, $query_pos, $query_seq, $matchlength, $upstream) =
    (0,'',0,0,'',0,0,0,0);
    
    my ($up_seq, $down_seq) = ('','');

    my $FH = $self->_output_fh;
    LINE:
    while (my $line = <$FH>) {
        chomp $line;

        if ( $line =~ m{MatchLength} ) {
            ( $matchlength, $subject_id, $subject_desc, $subject_pos, 
                $query_desc, $query_pos ) = $line =~ m{\A
                MatchLength:\s(\d+)\s
                "(.*?)\s(.*?)"
                \sat\s(\d+)\s
                vs\.\s
                "(.*?)"\s
                at\s(\d+)
             }xms;
             next LINE;   
        }    
        elsif ( $line =~ m{ \A > }xms ) { # -e mode
            ( $subject_id, $subject_desc, $subject_pos, 
                $query_desc, $query_pos ) = $line =~ m{\A
                >(.*?)\s
                (.*?)
                _at_(\d+)
                _with_
                (.*?)
                _at_(\d+)
             }xms;
             next LINE;   
        }    
        elsif ( $line =~ m{ \A 5 }xms) {
            ( $subject_seq ) = $line =~ m{ \A 5 (.*) 3 \z}xms;
            next LINE;
        }    
        elsif ( $line =~ m{ \A 3 }xms) {
            ( $query_seq ) = $line =~ m{ \A 3 (.*) 5 \z}xms;
            next LINE;
        }
        elsif ( $line =~ m{\A Maximum\snumber}xms ) {
            return 0;
        }    

        # find the query that belongs to the match
        my $query = $self->_mapping->{$query_desc};
        
        # already reverse, so just the complement here:
        $query_seq =~ tr[UTGCAutgca][AACGUaacgu];

        # -e mode
        if ( $line =~ m{ \A [gucaGUCA]+ \z }xms ) {
            # find subject in string. a little bit complicated because
            # don't know if the upstream/downstream region is as large
            # as we request
            my $qrc = lc($query->revcom->seq);
            $qrc =~ s/t/u/g;
            my $ql =  $s->query_length || length($s->query);
            SUBSTRING:
            for my $length ( reverse($ql .. length($s->query)) ) { 
                for my $start ( reverse( 0 .. $s->upstream )) {
                        my $query_start =
                        length($s->query)-$query_pos+1-$length;
                        my $qs = substr $qrc, $query_start, $length;
                        my $ss = substr $line, $start, $length;
                        #warn "L:$length S:$start QS:$query_start $qrc $qs $line $ss";
                        if ($self->_rnas_match($ss, $qs)) {
                                $subject_pos = $start;
                                $subject_seq = $ss;
                                $query_seq   = $qs;
                                $matchlength = $length;
                                $upstream    = $start;
                                $up_seq   = substr $line, 0, $upstream;
                                $down_seq = substr $line, $upstream + $length;
                                last SUBSTRING;
                        }        
                }        
            }        
        }    
        my $fasta = Bio::Seq->new(
            -id  => $subject_id,
            -seq => $up_seq .  $subject_seq . $down_seq,
            -desc => $subject_desc,
        );
        my $tmp_aln = new Bio::SimpleAlign( -source => "Bio::Grep" );
        $tmp_aln->add_seq(
            Bio::LocatableSeq->new(
                -id    => 'Subject',
                -seq   => $subject_seq,
                -start => $subject_pos,
                -end   => $subject_pos+ length($subject_seq)-1,
            )
        );
        my $rct = '';
        $rct =  '(Reverse Complement)' if $s->query_isset &&
            $s->reverse_complement;
        $tmp_aln->add_seq(
            Bio::LocatableSeq->new(
                -id    => "Query$rct",
                -seq   => $query_seq,
                -start => $query_pos,
                -end   => $query_pos + length($query_seq)-1,
            )
        );
        my $res = $self->_filter_result(
            Bio::Grep::Container::SearchResult->new( $fasta,
            $upstream, $upstream +$matchlength,
            $tmp_aln, $fasta->id, '' )
            );
        $res->query($query);    
        return $res if $res;
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

Bio::Grep::Backends::GUUGle - GUUGle back-end  


=head1 SYNOPSIS

  use Bio::Grep::Backends::GUUGle;
 
  use Bio::Root::Exception;
 
  # construct the GUUGle back-end	
  my $sbe = Bio::Grep::Backends::GUUGle->new();
 
  $sbe->settings->tmppath('tmp');
  $sbe->settings->datapath('data');
  
  # generate a GUUGle suffix array. you have to do this only once.
  $sbe->generate_database_out_of_fastafile('ATH1.cdna', 'AGI Transcripts (- introns, + UTRs)');
 
  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
 
  # take first available database in our test
  $sbe->settings->database($local_dbs[0]);
 
  my $seq = 'UGAACAGAAAGCUCAUGAGCC'; 
 
  # search for the reverse complement
  $sbe->settings->query($seq);
  $sbe->settings->reverse_complement(1);
  
  # display 5 bases upstream and downstream of the match
  $sbe->settings->upstream(5);
  # upstream and downstream must be the same in GUUGle
  $sbe->settings->downstream(5);
  
  $sbe->search();
 
  # output all informations we have!
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";
  }
   
=head1 DESCRIPTION

B<Bio::Grep::Backends::GUUGle> searches for a query in a GUUGle suffix array. 

WARNING: NOT THOROUGHLY TESTED! I recommend setting the environment variable BIOGREPDEBUG and
comparing the results with the program output.

NOTE 1: GUUGle always searches for the reverse complement. If you specify a
query file, Bio::Grep throws an exception if reverse_complement is not set.

NOTE 2: GUUGle only allows search for exact matches. It counts GU as no
mismatch.

NOTE 3: Upstream and Downstream must have the same length.


=head1 METHODS

See L<Bio::Grep::Backends::BackendI> for other methods. 

=over 2

=item Bio::Grep::Backends::GUUGle-E<gt>new()

    This function constructs a GUUGle back-end object

   my $sbe = Bio::Grep::Backends::GUUGle->new();

=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

   $sbe->sort('ga');

Available sortmodes in GUUGle:

=over

            ga  : 'ascending order of dG'
            gd  : 'descending order of dG'

=back

Note that 'ga' and 'gd' require that search results have dG set. 
L<Bio::Grep::RNA> ships with filters for free energy calculation. Also note that
these two sort options require that we load all results in memory.

=back


=head1 SEE ALSO

L<Bio::Grep::Backends::BackendI>
L<Bio::Grep::Container::SearchSettings>
L<Bio::Seq>


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
