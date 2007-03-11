package Bio::Grep::Backends::BackendI;

use strict;
use warnings;

use Bio::Grep::Container::SearchSettings;
use Bio::Grep::Root;

use Bio::AlignIO;
use Bio::Factory::EMBOSS;
use Bio::Index::Fasta;
use Bio::Seq;
use Bio::SeqIO;

use base 'Bio::Grep::Root';

use File::Spec;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;

use version; our $VERSION = qv('0.7.1');

use Class::MethodMaker [
    new      => 'new2',
    scalar   => [qw / settings _output_fh _output_fn _current_res_id/],
    array    => [qw / _query_seqs _results/],
    hash     => [qw / features _mapping/],
    abstract => [
        qw / search get_sequences get_databases
        generate_database_out_of_fastafile _parse_next_res/
    ],
];

sub new {
    my $self = shift->new2;

    # initialize standard settings
    my $settings = Bio::Grep::Container::SearchSettings->new();

    $self->settings($settings);

    # assume back-end binary is in path
    $self->settings->execpath('');

    my %all_features = (
        MISMATCHES        => 1,
        GUMISMATCHES      => 1,
        EDITDISTANCE      => 1,
        INSERTIONS        => 1,
        DELETIONS         => 1,
        FILTERS           => 1,
        NATIVE_ALIGNMENTS => 1,
        EVALUE            => 1,
        PERCENT_IDENTITY  => 1,
        PROTEINS          => 1,
        ONLINE            => 1,
        UPSTREAM          => 1,
        DOWNSTREAM        => 1,
        SORT              => 1,
        MAXHITS           => 1,
        COMPLETE          => 1,
        QUERYFILE         => 1,
        SHOWDESC          => 1,
        QSPEEDUP          => 1,
        REVCOM_DEFAULT    => 1,
    );
    $self->features(%all_features);
    $self;
}

sub _get_databases {
    my ( $self, $suffix ) = @_;
    my %result = ();
    opendir( DIR, $self->settings->datapath );
    my @files = grep {/${suffix}$/} readdir(DIR);
    closedir DIR;
    foreach my $file (@files) {
        my $prefix = $file;
        $prefix =~ s/$suffix//;
        $result{$prefix} = $prefix;
        $file =~ s/$suffix/\.nfo/;
        $file = $self->_cat_path_filename( $self->settings->datapath, $file );
        open my $FILE, '<', $file or next;
        my $desc = "";
        while ( my $line = <$FILE> ) {
            chomp $line;
            $desc .= $line;
        }
        $result{$prefix} = $desc;
        close $FILE
            or $self->throw(
            -class => 'Bio::Root::IOException',
            -text  => "Can't close $file" -value => $!
            );
    }
    return %result;
}

sub _filter_result {
    my ( $self, $res ) = @_;
    if ( $self->settings->filters_isset ) {
        foreach my $filter ( @{ $self->settings->filters } ) {
            $filter->message_reset;
            $res->_real_query( $self->settings->_real_query );
            $filter->search_result($res);
            unless ( $filter->filter ) {
                if ( $filter->delete ) {
                    return 0;
                }
            }
            $res->remark( $filter->message ) if $filter->message_isset;
        }
    }
    return $res;
}

sub results {
    my ( $self ) = @_;
    my @results;
    while (my $res = $self->next_res ) {
        push @results, $res;
    }
    wantarray ? @results : \@results;
}

sub _prepare_results {
    my ( $self ) = @_;
    if ( $self->settings->sort_isset
        && substr( $self->settings->sort, 0, 1 ) eq 'g' )
    {
        my @results = $self->results();
        @results = $self->_sort_by_dg(@results);
        $self->_results(@results);
    } else {
        $self->_results_reset;
    }
}

sub next_res {
    my ( $self ) = @_;
    if ( $self->_results_isset ) {
        my $id = $self->_current_res_id;
        $self->_current_res_id($id+1);
        return $self->_results->[$id];
    } else {
        return $self->_parse_next_res();
    }
}

sub _sort_by_dg {
    my ($self, @results) = @_;
    foreach my $result (@results) {
        if ( !$result->dG_isset ) {
            $self->warn(
                "Not sorting results by dG because some results have no
         dG calculated. Use Bio::Grep::RNA."
            );
            return @results;
        }
    }
    @results = sort { $a->dG <=> $b->dG } @results;
    if ( $self->settings->sort eq 'gd' ) {
        return reverse @results;
    }
    else {
        return @results;
    }
}

# calculates needleman-wunsch global alignment with the EMBOSS
# implementation
sub _get_alignment {
    my ( $self, $seq_a, $seq_b ) = @_;
    my $factory = Bio::Factory::EMBOSS->new();
    my $prog    = $factory->program('needle');
    my $outfile =
        $self->_cat_path_filename( $self->settings->tmppath,
        $seq_a->id . ".out" );
    my @seqs    = qw($seq_b);
    my $gapopen = '5.0';
    $gapopen = '50.0' unless $self->settings->editdistance_isset;

    $prog->run(
        {   -asequence => $seq_a,
            -bsequence => $seq_b,
            -gapopen   => $gapopen,
            -gapextend => '5.0',
            -outfile   => $outfile
        }
    );

    if ( -e $outfile ) {
        my $alignio_fmt = "emboss";
        my $align_io    = new Bio::AlignIO(
            -format => $alignio_fmt,
            -file   => $outfile
        );
        unlink($outfile);
        return $align_io->next_aln();
    }

}

# concatenates a path and a filename platform-independently
sub _cat_path_filename {
    my ( $self, $path, $filename ) = @_;
    return $filename if $path eq '';
    return File::Spec->catfile( $path, $filename );
}

sub _check_search_settings {
    my ($self, $arg_ref) = @_;
    
    if (defined $arg_ref) {
        $self->settings->set($arg_ref);
    }    
    $self->settings->upstream(0)   if !defined $self->settings->upstream;
    $self->settings->downstream(0) if !defined $self->settings->downstream;

    $self->settings->upstream(
        $self->is_integer( $self->settings->upstream ) );
    $self->settings->downstream(
        $self->is_integer( $self->settings->downstream ) );

    $self->settings->mismatches(
        $self->is_integer( $self->settings->mismatches ) );
    $self->settings->mismatches_reset if !defined $self->settings->mismatches;

    $self->settings->insertions(
        $self->is_integer( $self->settings->insertions ) );
    $self->settings->insertions_reset if !defined $self->settings->insertions;

    $self->settings->deletions(
        $self->is_integer( $self->settings->deletions ) );
    $self->settings->deletions_reset if !defined $self->settings->deletions;
    
    $self->settings->showdesc(
        $self->is_integer( $self->settings->showdesc ) );
    $self->settings->showdesc_reset if !defined $self->settings->showdesc;
    
    $self->settings->qspeedup(
        $self->is_integer( $self->settings->qspeedup ) );
    $self->settings->qspeedup_reset if !defined $self->settings->qspeedup;

    $self->settings->editdistance(
        $self->is_integer( $self->settings->editdistance ) );
    $self->settings->editdistance_reset
        if !defined $self->settings->editdistance;

    $self->settings->query_length(
        $self->is_integer( $self->settings->query_length ) );
    $self->settings->query_length_reset
        if !defined $self->settings->query_length;

    $self->settings->maxhits( $self->is_integer( $self->settings->maxhits ) );
    $self->settings->maxhits_reset if !defined $self->settings->maxhits;

    if ( $self->settings->sort_isset ) {
        my $found_sort_mode = 0;
        my %sort_modes      = $self->available_sort_modes();
        foreach my $sort_mode ( keys %sort_modes ) {
            if ( $self->settings->sort eq $sort_mode ) {
                $sort_mode =~ /(\w+)/;
                $self->settings->sort($1);    #make taint happy
                $found_sort_mode = 1;
                last;
            }
        }
        if ( $found_sort_mode == 0 ) {
            $self->throw(
                -class => 'Bio::Root::BadParameter',
                -text  => 'Sort mode not valid',
                -value => 'sort mode'
            );
        }

    }

    # check if database is set and valid
    my $found_database = 0;
    if ( defined $self->settings->database ) {
        $self->settings->database(
            $self->is_word( $self->settings->database ) );
    }
    else {
        $self->throw(
            -class => 'Bio::Root::BadParameter',
            -text  => 'Database not defined',
        );
    }

    # some warnings if user requests features that are not available
    unless ( defined( $self->features->{GUMISMATCHES} ) ) {
        $self->warn("GU mismatch options not available in this back-end")
            if $self->settings->gumismatches_isset
            && $self->settings->gumismatches != 1;
    }
    unless ( defined( $self->features->{DELETIONS} ) ) {
        $self->warn(
            "deletions not available in this back-end. Try editdistance or mismatches."
            )
            if $self->settings->deletions_isset
            && $self->settings->deletions != 0;
    }
    unless ( defined( $self->features->{INSERTIONS} ) ) {
        $self->warn(
            "insertions not available in this back-end. Try editdistance or mismatches."
            )
            if $self->settings->insertions_isset
            && $self->settings->insertions != 0;
    }
    unless ( defined( $self->features->{FILTERS} ) ) {
        $self->warn("Filter options not available in this back-end")
            if $self->settings->filters_isset;
    }
    unless ( defined( $self->features->{ONLINE} ) ) {
        $self->warn("online parameter not available in this back-end")
            if $self->settings->online_isset && $self->settings->online;
    }
    unless ( defined( $self->features->{COMPLETE} ) ) {
        $self->warn("complete parameter not available in this back-end")
            if $self->settings->complete_isset && $self->settings->complete;
    }
    unless ( defined( $self->features->{SHOWDESC} ) ) {
        $self->warn("showdesc parameter not available in this back-end")
            if $self->settings->showdesc_isset;
    }
    unless ( defined( $self->features->{QSPEEDUP} ) ) {
        $self->warn("qspeedup parameter not available in this back-end")
            if $self->settings->qspeedup;
    }
    if (defined $self->settings->editdistance && defined $self->settings->mismatches &&
    $self->settings->editdistance > 0 && $self->settings->mismatches > 0) {
        $self->settings->editdistance_reset;
        $self->warn(
            "Ignoring edit distance settings because mismatches is set.")
          ;
    }        
    if ( $ENV{BIOGREPDEBUG} ) {
        warn $self->settings->to_string();
    }

    if ( $self->settings->filters_isset ) {
        foreach my $filter ( @{ $self->settings->filters } ) {
            $filter->reset;
        }
    }
    $self->_results_reset;
    $self->_current_res_id(0);
    return;
}

# copies the specified fasta file in the data directory and cerate a file
# <databasename>.nfo with the description, specified in the optional 2nd
# argument
# changes directory! so please save the oldpath before calling this function
sub _copy_fasta_file_and_create_nfo {
    my ( $self, $file, $filename, $description ) = @_;
    
    # throw exception if filename looks wrong
    $self->is_word($filename, 'FASTA filename');
    
    my $newfile =
        $self->_cat_path_filename( $self->settings->datapath, $filename );

    copy( $file, $newfile )
        or $self->throw(
        -class => 'Bio::Root::IOException',
        -text  => "Can't copy $file to $newfile",
        -value => $!
        );

    if ( defined($description) ) {
        open my $NFOFILE, '>', $newfile
            . '.nfo'
            or $self->throw(
            -class => 'Bio::Root::FileOpenException',
            -text  => "Can't open $filename.nfo for writing" -value => $!
            );
        print $NFOFILE $description
            or $self->throw(
            -class => 'Bio::Root::IOException',
            -text  => "Can't write to $filename.nfo" -value => $!
            );
        close $NFOFILE
            or $self->throw(
            -class => 'Bio::Root::IOException',
            -text  => "Can't close $filename.nfo" -value => $!
            );

    }
}

sub _guess_alphabet_of_file {
    my ( $self, $filename ) = @_;
    my $in = Bio::SeqIO->new( -file => $filename, -format => 'fasta');
    return $in->next_seq->alphabet;
}

# prepares the query, for example calculating the reverse complement if necessary
# returns the prepared query. settings->query is unchanged!
sub _prepare_query {
    my $self  = shift;
    my $query = $self->settings->query;
    my $db_alphabet =
    $self->get_alphabet_of_database($self->settings->database);
    
    if ($db_alphabet eq 'dna') {
        $query =~ tr/uU/tT/;
    }
    my $seq = Bio::Seq->new( -seq => $query );
    if ($seq->alphabet ne $db_alphabet) {
        $self->throw( -class => 'Bio::Root::BadParameter',
                      -text  => 'Alphabet of query and database not equal',
                      -value => 'Seq: '. $seq->alphabet . ", DB: $db_alphabet"
                    );
    }
    $self->settings->reverse_complement(0)
        unless $self->settings->reverse_complement_isset;
    if ($self->settings->reverse_complement) {
        if ($db_alphabet eq 'dna') {
            $query = $seq->revcom->seq unless defined 
                $self->features->{REVCOM_DEFAULT}; 
        }
        else {
            $self->warn("Setting reverse complement only available for DNA
            databases. I will unset reverse complement");
            $self->settings->reverse_complement(0);
        }
    }
    elsif (defined $self->features->{REVCOM_DEFAULT} && $db_alphabet eq 'dna') {    
            $query = $seq->revcom->seq 
    }    
    $self->settings->_real_query( uc($query) );
    return $self->settings->_real_query();
}

sub available_sort_modes {
    my ($self) = @_;
    return (
        ga => 'ascending order of dG',
        gd => 'descending order of dG',
    );
}

sub _get_sequences_from_bio_index {
    my ( $self, $seqid ) = @_;
    my $indexfile =
        $self->settings->datapath . '/' . $self->settings->database . '.idx';
    my $idx = Bio::Index::Fasta->new($indexfile);
    my $string;
    my $stringio = IO::String->new($string);
    my $out = Bio::SeqIO->new( -fh => $stringio, -format => 'Fasta' );
    foreach my $seqid ( @{$seqid} ) {
        $out->write_seq( $idx->fetch($seqid) );
    }
    $stringio = IO::String->new($string);
    return Bio::SeqIO->new( -fh => $stringio, -format => 'Fasta' );
}

sub get_alphabet_of_database {
    my ( $self, $db ) = @_;
    my %dbs = $self->get_databases();
    if (!defined $dbs{$db}) {
        my $self->throw( -class => 'Bio::Root::BadParameter',
                         -text  => 'Database not found',
                         -value => $db,
                       );
    }
    $db = $self->is_word($db);
    
    open my $ALFILE, '<', $self->_cat_path_filename($self->settings->datapath,
    $db . '.al1');
    my $lines = 0;
    while (my $line = <$ALFILE>) {
        $lines++;
    }
    close $ALFILE;
    $lines <= 5 ? return 'dna' : return 'protein';
}

sub _delete_output {
    my ( $self ) = @_;
    return 0 if !defined $self->_output_fn;
    return 0 if !-e $self->_output_fn;
    unlink $self->_output_fn or
        $self->throw(
            -class => 'Bio::Root::IOException',
            -text  =>
                "Cannot remove " . $self->_output_fn,
            -value => $!,    
            );
   return 1;
}

sub _execute_command {
    my ( $self, $cmd ) = @_;
    $self->_delete_output();
    my ( $tmp_fh, $tmp_fn ) =
            tempfile( 'parser_XXXXXXXXXXXX', DIR => $self->settings->tmppath );
    $self->_output_fh($tmp_fh);        
    $self->_output_fn($tmp_fn);        
    system("$cmd > $tmp_fn");
    
    return !$?;
}

sub _create_index_and_alphabet_file {
    my ( $self, $filename ) = @_;
    my $in = Bio::SeqIO->new( -file => $filename, -format => 'Fasta' );
    my $alphabet = $in->next_seq()->alphabet();
    if (!defined $self->features->{PROTEINS} && $alphabet eq 'protein') {
        $self->throw(
        -class => 'Bio::Root::BadParameter',
        -text  => "Back-end does not support protein data",
        -value => $alphabet,
        );
    }    
    # create a vmatch alphabet file 
    open my $ALFILE, '>', "$filename.al1"
        or $self->throw(
        -class => 'Bio::Root::FileOpenException',
        -text  => "Can't open $filename.al1 for writing",
        -value => $!
        );
    if ($alphabet eq 'dna') {
        print $ALFILE "aA\ncC\ngG\ntTuU\nnsywrkvbdhmNSYWRKVBDHM\n";
    }
    else {
        print $ALFILE
        "L\nV\nI\nF\nK\nR\nE\nD\nA\nG\nS\nT\nN\nQ\nY\nW\nP\nH\nM\nC\nXUBZ*\n";
    }
    close $ALFILE;
    # create database from directory of fasta files
    my $idx = Bio::Index::Fasta->new(
        -filename   => $filename . '.idx',
        -write_flag => 1
    );
    $idx->make_index( ($filename) );

    return;
}

sub _create_tmp_query_file {
    my ( $self ) = @_;
    my $s = $self->settings;
    my $query_file = 0;
    my ( $tmp_fh, $tmp_query_file );
    if (defined $s->query_file) {
        $query_file = $self->is_path($s->query_file);
    }
    my @query_seqs = ();
    
    my $query = '';
    if (!$query_file) {
        $query = $self->_prepare_query();

        ( $tmp_fh, $tmp_query_file ) =
            tempfile( 'spatter_XXXXXXXXXXXX', DIR => $s->tmppath, UNLINK =>
                !$ENV{BIOGREPDEBUG} ); # don't delete when in debug mode, so
                                       # we can reproduce it.

        # construct a temporary fasta file with the query for vmatch
        my $seqobj = Bio::Seq->new(
            -id => 'Query',
            -desc       => $tmp_query_file,
            -seq        => $query
        );

        my $outseqio = Bio::SeqIO->new(
            -fh     => $tmp_fh,
            -format => 'fasta'
        );
        $outseqio->write_seq($seqobj);
        push @query_seqs, $seqobj;
    }
    else {
        $tmp_query_file = $query_file;
        my $query_in = Bio::SeqIO->new(-file => $tmp_query_file);
        while ( my $query_seq = $query_in->next_seq) {
            push @query_seqs, $query_seq;
        }
    }
    $self->_query_seqs(@query_seqs);
    return ( $query, $query_file, $tmp_query_file);
}   

1;    # Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep::Backends::BackendI - Superclass for all back-ends  

=head1 SYNOPSIS

See the back-end modules for example code.

=head1 DESCRIPTION

B<Bio::Grep::Backends::BackendI> is the superclass for all back-ends. Don't use this class
directly. 

=head1 METHODS

=over 

=item C<new()>

This function constructs a Backend object.

=cut

=item C<$sbe-E<gt>settings()>

Get the settings. This is a L<Bio::Grep::Container::SearchSettings> object

  # search for the reverse complement and allow 4 mismatches
  $sbe->settings->database('ATH1.cdna');
  $sbe->settings->query('UGAACAGAAAGCUCAUGAGCC');
  $sbe->settings->reverse_complement(1);
  $sbe->settings->mismatches(4);

=item C<$sbe-E<gt>results()>

Get the results. This is an array of L<Bio::Grep::Container::SearchResults> objects.

  # output the searchresults with alignments
  foreach my $res (@{$sbe->results}) {
     print $res->sequence->id . "\n";
     print $res->alignment_string() . "\n\n";
  }

This method is DEPRECATED. The new syntax is

  # output the searchresults with alignments
  while ( my $res = $sbe->next_res ) {
     print $res->sequence->id . "\n";
     print $res->alignment_string() . "\n\n";
  }
 
  if you need an array with all search results, you should use following code:

  my @results;

  while ( my $res = $sbe->next_res ) {
      push @results, $res;
  }
          

=item C<$sbe-E<gt>features()>

Get available features. This is a hash. Valid features are
MISMATCHES, GUMISMATCHES, EDITDISTANCE, INSERTIONS, DELETIONS, 
FILTERS, NATIVE_ALIGNMENTS, PROTEINS, UPSTREAM, DOWNSTREAM, MAXHITS, COMPLETE,
QUERYFILE, SHOWDESC, QSPEEDUP, EVALUE and PERCENT_IDENTITY.

  if (defined($sbe->features->{GUMISMATCHES})) {
          # $sbe->settings->gumismatches(0);
          $sbe->settings->gumismatches(0.5);
  } else {
        print "\nBack-end does not support wobble pairs\n";
  }

=item C<$sbe-E<gt>get_alphabet_of_database($db)>

Returns 'dna' if the specified database is a DNA database, 'protein'
otherwise.

=back

=head1 ABSTRACT METHODS

Every back-end must implement this methods.

=over

=item C<$sbe-E<gt>search>

This function searches for the query specified in the 
L<Bio::Grep::Container::SearchSettings> object 
C<$sbe-E<gt>settings>
 
    $sbe->search();

=item C<$sbe-E<gt>next_res>

Returns next result. L<Bio::Grep::Container::SearchResult> object

    while ( my $res = $sbe->next_res ) {
        # output result
    }
    
=item C<$sbe-E<gt>get_sequences>

This function returns all sequences with the ids in the specified array reference as
a L<Bio::SeqIO> object. 


   my $seqio = $sbe->get_sequences([$id]);
   my $string;  my $stringio = IO::String->new($string);
   my $out = Bio::SeqIO->new('-fh' => $stringio,
                             '-format' => 'fasta');

   while ( my $seq = $seqio->next_seq() ) {
      # write the sequences in a string
      $out->write_seq($seq);
   }
   print $string;
 

=item C<$sbe-E<gt>get_databases>

Returns a hash with all available databases. The keys are the filenames,
the values are descriptions (or the filename if no description is available).

Descriptions can be set in info files. For example, if you indexed file
ATH1.cdna, Vmatch and HyPA construct a lot of ATH1.cdna.* files. Now simply create a file 
ATH1.cdna.nfo and write a description in that file. The function C<generate_database_out_of_fastafile> will
create this file for you if you add a description as second argument.

  my %local_dbs_description = $sbe->get_databases();
  my @local_dbs = sort keys %local_dbs_description;
  
  # take first available database 
  $sbe->settings->database($local_dbs[0]);


=item C<$sbe-E<gt>generate_database_out_of_fastafile($fastafile)>

Copies the specified file in the datapath directory (C<$sbe-E<gt>settings-E<gt>datapath>)
and generates a database (HyPa/Vmatch: a suffix array). You can get the
available databases with C<$sbe-E<gt>get_databases()>. You have to do this
only once. Vmatch and HyPa need a lot of RAM for the construction of their enhanced
suffix arrays.

  $sbe->generate_database_out_of_fastafile('ATH1.cdna', 'AGI Transcripts');

=item C<$sbe-E<gt>available_sort_modes()>

Returns a hash with the available result sort modes. Keys are the modes you
can set with $sbe->settings->sort($mode), values a short description.

=back

=head1 INTERNAL METHODS 

Only back-ends should call them directly.

=over 

=item C<_check_search_settings> 

Performs some basic error checking. Important security checks, because
we use system(). So we should check, if we get what we assume.

Because every back-end should call this function at the top 
of its search function, we clean things like old search results here up

=item C<_prepare_query>

Another important method that every back-end must call.
Prepares the query, for example calculating the reverse complement if
necessary, returns the prepared query. settings->query is unchanged!


=item C<_copy_fasta_file_and_create_nfo>

The generate_database_out_of_fastafile implementation of your back-end class
should use this function to copy the specified Fasta file to the data
directory  and to generate an info file, containing the description of the
Fasta file.

=item C<_get_alignment( $seq_query, $seq_subject )>

Calculates and returns an alignment of two L<Bio::Seq> objects. Requires
EMBOSS and bioperl-run.

=item C<_get_databases($suffix)>

This function searches the data directory for files ending with $suffix
and returns this list of files in an array

Substitutes $suffix with .nfo from all found files and searches for an
info file with that name. The content of that file will be used as description.
When no file is found, the description will be the filename without the suffix:


   %dbs = _get_databases('.al1');  # finds file ATH1.cdna.al1, searches for ATH1.cdna.nfo
   print $dbs{'ATH1.cdna'};      # prints content of ATH1.cdna.nfo or 'ATH1.cdna'


=item C<_get_sequences_from_bio_index($id)>

Hypa and Agrep back-end use L<Bio::Index> for sequence id queries (implemented
in this this method. Returns a L<Bio::SeqIO> object like abstract the method get_sequences should.

=item C<_create_tmp_query_file()>

Examines query, query_file and reverse_complement and generates a temporary
Fasta file that is passed in the C<system()> call to the backend. If the
environment variable BIOGREPDEBUG is not set, then this file will be deleted
when the script exits.  

=item C<_create_index_and_alphabet_file($fastafile)>

Creates an index of the specified Fasta file with L<Bio::Index::Fasta>.
Creates an Vmatch alphabet file.

=back

=head1 DIAGNOSTICS

=over

=item C<Bio::Root::IOException>

It was not possible to copy the Fasta file in
C<generate_database_out_of_fastafile> into the data directory or
it was not possible to write to the Fasta file in the data directory.
Check permissions, data path and free disk space.

=item C<Bio::Root::FileOpenException>

It was not possible to open the Fasta file in the data directory for
writing. Check permissions.

=item C<Bio::Root::BadParameter>

You started the search with invalid search settings. 

=over

=item C<Sort mode not valid>

The specified sort mode ($sbe->settings->sort) is not valid.
You can get all valid sort modes with $sbe->available_sort_modes()
See L<Bio::Grep::Backends::Vmatch>, L<Bio::Grep::Backends::Hypa>,
L<Bio::Grep::Backends::Agrep> for details.

=item C<Database not defined>

You forgot to define a database. You have to build a database with
$sbe->generate_database_out_of_fastafile (once) and set it with
$sbe->settings->database(). Example:

$sbe->generate_database_out_of_fastafile('ATH1.cdna");
$sbe->settings->database('ATH1.cdna');
   


=item C<Database not valid (insecure characters)>

The database name is not valid. Allowed characters are 'a-z', 'A-z','0-9', '.'
, '-' and '_'.


=back

=back


=head1 FILES

Requires EMBOSS and Bio::Factory::EMBOSS for the Needleman-Wunsch local 
alignment implementation from EMBOSS. The internal method 
C<_get_alignment($seq_a, $seq_b)> can than calculate an alignment for 
back-ends that do not generate a alignment (like Hypa, agrep).


=head1 SEE ALSO

L<Bio::Grep::Container::SearchSettings> 
L<Bio::Grep::Container::SearchResults> 


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
