package Bio::Grep::Backends::Vmatch;

use strict;
use warnings;

use Bio::Grep::Container::SearchResult;
use Bio::Grep::Backends::BackendI;

use base 'Bio::Grep::Backends::BackendI';

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use IO::String;

use version; our $VERSION = qv('0.9.3');

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
    my ($self, $arg_ref) = @_;
    my $s    = $self->settings;
    $self->_check_search_settings($arg_ref);

    my ( $query, $query_file, $tmp_query_file) =
    $self->_create_tmp_query_file();

    if (($s->upstream > 0 || $s->downstream > 0) && $s->showdesc_isset) {
        $self->throw(-class => 'Bio::Root::BadParameter',
              -text  => "You can't use showdesc() with upstream or downstream.",
                );
    }
    if ($query_file && !$s->complete && !$s->query_length_isset) {
        $self->throw(-class => 'Bio::Root::BadParameter',
              -text  => 'You have to specify complete or querylength. See ' .
             'the flags -complete and -l in the Vmatch documentation.',
                );
    }    
    # now generate the command string
    my $fuzzy = '';
    $fuzzy = ' -h ' . $s->mismatches . ' ' if $s->mismatches > 0;
    if ( $s->editdistance_isset ) {
        $fuzzy = ' -e ' . $s->editdistance . ' ';
    }
    my $online = '';
    $online = ' -online ' if ( $s->online_isset && $s->online );
    my $auto_query_length = 0;
    if (!defined $s->query_length && !$s->complete && !$query_file) {
        $s->query_length( length($query) );
        $auto_query_length = 1;
    }    
    my $length = '';
    $length =  ' -l ' . $s->query_length if defined $s->query_length;
    my $complete = '';
    $complete = ' -complete ' if ($s->complete_isset && $s->complete);
    
    my $sort   = '';

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
            -text  =>
                "Vmatch: You can't combine qspeedup and complete",
            ) if $s->complete_isset;

        $qspeedup = ' -qspeedup ' . $s->qspeedup . ' '; 

    }

    my $hxdrop = '';
    $hxdrop = ' -hxdrop ' . $s->hxdrop . ' ' if $s->hxdrop_isset;
    
    my $exdrop = '';
    $exdrop = ' -exdrop ' . $s->exdrop . ' ' if $s->exdrop_isset;
   
    my $pflag = '';
    $pflag = ' -p ' if ($s->query_file && $s->reverse_complement);

    my $command =
        $self->_cat_path_filename( $s->execpath, 'vmatch' ) . ' -q '
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
        -text  =>
            "Vmatch error: Query not valid(maybe too many mismatches? Try online(1))\n\nCommand was:\n\t$command"
        )
        if !$cmd_ok;
        
#    warn $output;   
    $self->_prepare_results;
    
    if ($s->showdesc_isset) {
        my %query_desc_lookup;
        foreach my $query_seq (@{ $self->_query_seqs }) {
            # simulate how this sequence would look in vmatch output
            my $query_desc = $query_seq->id;
            $query_desc .= '_' . $query_seq->desc if length $query_seq->desc;
            $query_desc = substr $query_desc, 0, $s->showdesc;
            $query_desc =~ s{ }{_}g;
            $query_desc_lookup{$query_desc} = $query_seq;
        }
        $self->_mapping(%query_desc_lookup);
        #warn Data::Dumper->Dump([%query_desc_lookup]);
    }
    $self->settings->query_length_reset if $auto_query_length;
    return 1;
}

sub get_databases {
    my $self = shift;
    return $self->_get_databases('.al1');
}


sub generate_database_out_of_fastafile {
    my $self        = shift;
    my $file        = shift;
    my $description = shift;
    my ( $filename, $oldpath ) = fileparse($file);

#    my $olddir = getcwd();

    $self->_copy_fasta_file_and_create_nfo( $file, $filename, $description );
    
    my $alphabet = $self->_guess_alphabet_of_file($file);
    my $alphabet_specific_arguments = '';
    #warn $alphabet; 
    if ($alphabet eq 'protein') {
        $alphabet_specific_arguments = ' -protein -pl -allout -v ';
    }
    elsif ($alphabet eq 'dna') {
        $alphabet_specific_arguments = ' -dna -pl 3 -allout -v ';
    }
    else {
        $self->throw(-class => 'Bio::Root::BadParameter',
                     -text  => 'unsupported alphabet of file',
                     -value => $alphabet,);   
    }
    my $command =
        $self->_cat_path_filename( $self->settings->execpath, 'mkvtree' )
        . ' -db '
        . $filename
        . $alphabet_specific_arguments;
        
    #warn $command;
    my $output_dir = $self->settings->datapath;
    system(qq{ cd $output_dir ; exec $command } );
    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  =>
            "Vmatch error: Cannot generate suffix array. Command was:\n\t$command"
        )
        if ($?);
   return $filename;     
}

sub _parse_next_res {
    my $self                = shift;
    my @query_seqs = $self->_query_seqs;
    my $s                   = $self->settings;
    my @results             = ();
    my $alignment_in_output = -1;
    my $skip_next_alignment = 0;
    my $skipped_lines = 0;
    my $subject;
    my $tmp_aln;
    
    my ($command, $output);  

    my $FH = $self->_output_fh;
    LINE:
    while (my $line = <$FH>) { 
        chomp $line;
        $line =~ s/\s+/ /g;
        if ($line !~ /\s/) {
            $skipped_lines++;
            if ($skipped_lines == 2) {
                $results[-1]->alignment($tmp_aln);
                my $real_subject = $subject;
                # remove gaps out of alignment
                $real_subject =~ s{-}{}g;
                $results[-1]->sequence->seq($real_subject) if $s->showdesc_isset;
    #            $results[-1]->query->seq($query);
                my $res = $self->_filter_result($results[-1]);
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
        $alignment_in_output = 0 if ( $line =~ /^Sbjct:/ && !$skip_next_alignment );

        next unless ( $fields[0] =~ /^\d+$/ || $alignment_in_output >= 0 );
        $skip_next_alignment = 0;

        if ($line =~ /^Sbjct: (.*)$/) {
            $subject = $1;
        }
        if ($line =~ /^Query: (.*)$/) {
            my $query = $1;
            $query =~ s/\s+(\d+)\s*$//;
            my $query_pos = $1;
            $subject =~ s/\s+(\d+)\s*$//;
            my $subject_pos = $1;
            if (!$tmp_aln->no_sequences) {
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
                $s1->seq($s1->seq . $subject);
                $s2->end($query_pos);
                $s2->seq($s2->seq . $query);
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

        if ($s->showdesc_isset) {
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
        my $upstream = $s->upstream;
        my $internal_seq_id = $fields[1];
        
        if (!$s->showdesc_isset) {        
            my $start = $fields[2] - $s->upstream;

            # maybe the defined upstream region is larger than available
            # so check this and store in local variables
            if ( $start < 0 ) {
                $upstream = $upstream + $start;
                $start    = 0;
            }
            my $length  = $upstream + $fields[0] + $s->downstream;
            $command =
                $self->_cat_path_filename( $s->execpath, 'vsubseqselect' )
                . " -seq $length $fields[1] "
                . $start . ' '
                . $self->_cat_path_filename( $s->datapath, $s->database );
            $output   = `$command`;
            my $stringio = IO::String->new($output);
            my $in       = Bio::SeqIO->new(
                '-fh'     => $stringio,
                '-format' => 'fasta'
            );
            $fasta = $in->next_seq();
        }
        else {
            my ($seq_id, $seq_desc) = $fields[1] =~ m{\A (.+?) _ (.*) \z}xms;
            $seq_id = $fields[1] if !defined $seq_id;
            $seq_desc = '' if !defined $seq_desc;
            $seq_desc =~ s{_}{ }g;
            $fasta = Bio::Seq->new(-id => $seq_id, -desc => $seq_desc);
            $internal_seq_id = $seq_id;
        }
#        warn Data::Dumper->Dump([@lines]) if !defined $fasta->seq;
        
        my $alignment =
            Bio::SimpleAlign->new()
            ;    #$self->_get_alignment($seq_query, $seq_subject);
        my $result =
            Bio::Grep::Container::SearchResult->new( $fasta, $upstream,
            $upstream + $fields[0],
            $alignment, $internal_seq_id, '' );
        $result->evalue( $fields[8] );
        $result->percent_identity( $fields[10] );
        if ($s->showdesc_isset) {
            $result->query( $self->_mapping->{$fields[5]} );
        } 
        else {    
            $result->query( $query_seqs[$fields[5]] );
        }
        push ( @results, $result );
    }
    $self->_delete_output();
    return 0;
}


sub get_sequences {
    my ( $self, $seqid)  = @_;
    my $s     = $self->settings;
    $self->_check_search_settings();
    my ( $tmp_fh, $tmpfile );
    
    my $seq_query = '';
    
    if (@{$seqid}[0] =~ m{\A \d+ \z}xms) {
        ( $tmp_fh, $tmpfile ) =
            tempfile( 'vseqselect_XXXXXXXXXXXXX', DIR => $s->tmppath );

        foreach ( @{$seqid} ) {
            print $tmp_fh $_ . " \n ";
        }
        close $tmp_fh;
        $seq_query = ' -seqnum ' . $tmpfile;    
    }
    else {
        my $seq_desc = $self->is_sentence(@{$seqid}[0]);
#        ( $seq_desc) = $seq_desc =~ m{ \A (.*?) _? (.*) \z }xms;        
        $seq_query = ' -matchdesc "' . $seq_desc . '"';
    }

    my $command =
        $self->_cat_path_filename( $s->execpath, 'vseqselect' )
        . $seq_query . ' ' 
        . $self->_cat_path_filename( $s->datapath, $s->database );
#warn $command;
    my $output = `$command`;
    if ($? && $output !~ m{\A \> }xms) {
    $self->throw(
        -class => 'Bio::Root::SystemException',
        -text  =>
            "Vmatch error: Cannot fetch sequence out of suffix array. Command was:\n\t$command\n$output"
        )
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

Bio::Grep::Backends::Vmatch - Vmatch back-end  

=head1 SYNOPSIS

  use Bio::Grep::Backends::Vmatch;
  
  use Bio::Root::Exception;
  use Error qw(:try);
  
  # configure our search back-end, in this case Vmatch
  my $sbe = Bio::Grep::Backends::Vmatch->new();
  
  $sbe->settings->execpath('/usr/local/vmatch');
  $sbe->settings->tmppath('/tmp');
  $sbe->settings->datapath('data');
  
  # generate a Vmatch suffix array. you have to do this only once.
  $sbe->generate_database_out_of_fastafile('ATH1.cdna', 'AGI Transcripts (- introns, + UTRs)');
  
  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
  
  # take first available database in our test
  $sbe->settings->database($local_dbs[0]);
  
  # search for the reverse complement and allow 4 mismatches
  $sbe->settings->query('UGAACAGAAAGCUCAUGAGCC');
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(4);
  
  # With many mismatches and short queries, the "online" algorithm
  # is maybe faster. This alogrithm does not use the index. Test this!
  # $sbe->settings->online(1);

  # if you don't need upstream/downstream options, set showdesc
  # (see vmatch manual) for performance reasons 
  $sbe->settings->showdesc(100)
  
  $sbe->search();

  # output the searchresults with nice alignments
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->mark_subject_uppercase() . "\n";
     print $res->alignment_string() . "\n\n";
  }
  
=head1 DESCRIPTION

B<Bio::Grep::Backends::Vmatch> searches for a query in a Vmatch suffix array. 

=head1 METHODS

See L<Bio::Grep::Backends::BackendI> for other methods. 

=over 2

=item Bio::Grep::Backends::Vmatch-E<gt>new()

This function constructs a Vmatch back-end object

   my $sbe = Bio::Grep::Backends::Vmatch->new();


=item C<$sbe-E<gt>available_sort_modes()>

Returns all available sort modes as hash. keys are sort modes, values a short
description.

   $sbe->sort('ga');

Available sortmodes in Vmatch:

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
L<Bio::Grep::RNA> ships with filters for free energy calculation.  Also note that
these two sort options require that we load all results in memory.

=item C<$sbe-E<gt>get_sequences()>

Takes as argument an array reference. If first array element is an integer, then
this method assumes that the specified sequence ids are Vmatch internal ids.
Otherwise it will take the first array element as query.

    # get sequences 0,2 and 4 out of suffix array
    $sbe->get_sequences([0,2,4]);

    # get sequences that start with At1g1
    $sbe->get_sequences(['At1g1', 'ignored']);


=back

=head1 DIAGNOSTICS

See L<Bio::Grep::Backends::BackendI> for other diagnostics. 

=over

=item C<Bio::Root::SystemException>

=over 2

=item C<Vmatch error: Query not valid ... > 

It was not possible to run Vmatch in function C<search>. Check the search
settings.

=item C<Vmatch error: Cannot generate suffix array> 

It was not possible to generate a suffix array in function
C<generate_database_out_of_fastafile>. Check permissions and paths.

=item C<Vmatch error: Cannot fetch sequence out of suffix array> 

It was not possible to get some sequences out of suffix array in function
C<get_sequences>. Check sequence ids.

=back

=item C<Bio::Root::BadParameter>

=over 2

=item C<You can't use showdesc() with upstream or downstream.>

We need the tool C<vsubseqselect> of the Vmatch package for the upstream and 
downstream regions. That tools requires as parameter an internal vmatch
sequence id, which is not shown in the Vmatch output when showdesc is on.

=item C<You have to specify complete or querylength. ...'>

The Vmatch parameters -complete and -l cannot combined. See the Vmatch 
documentation.

=item C<Vmatch: You can't combine qspeedup and complete>

The Vmatch parameters -complete and -qspeedup cannot combined. See the Vmatch 
documentation.

=item C<unsupported alphabet of file>

The method generate_database_out_of_fastafile() could not determine the
alphabet (DNA or Protein) of the specified Fasta file.

=back

=item It is too slow, if I call vmatch on the command line, it is much faster!

Did you set showdesc(100)? Yes? Write a bugreport!

=back

=head1 SEE ALSO

L<Bio::Grep::Backends::BackendI>
L<Bio::Grep::Container::SearchSettings>
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
