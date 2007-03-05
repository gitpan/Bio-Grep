package BioGrepTest;

use File::Spec;


sub find_binary_in_path {
   my ( $name ) = @_;
   my @PATH = File::Spec->path();
   foreach my $path (@PATH) {
      if ( -e $path . "/$name" && !-d $path . "/$name") {
         return $path ;
      }
   }
   return '';
}

sub set_path {
    my ( @execs) = @_;
    my @paths = ( '/bin', '/usr/bin' );
    
    foreach my $exec (@execs) {
        my $path = find_binary_in_path($exec);
        push @paths, $path if $path ne '';
    }
    ( $ENV{PATH} ) = join(':', @paths) =~ /(.*)/; 
}

sub delete_files {
    # delete everything in this directory, otherwise we fail the test
    foreach my $file ( <t/tmp/*> ) {
        my ( $filename ) = $file =~ m{\A t/tmp/ ([\d\w\.\-]+) \z}xms;
        warn $file if !defined $filename;
        unlink "t/tmp/$filename";
    }
    foreach my $file ( <t/data/*> ) {
        my ( $filename ) = $file =~ m{\A t/data/ ([\d\w\.\-]+) \z}xms;
        warn $file if !defined $filename;
        unlink "t/data/$filename";
    }
    return 1;
}        

1;
