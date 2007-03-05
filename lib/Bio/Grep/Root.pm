package Bio::Grep::Root;

use strict;
use warnings;

use Bio::Root::Root;

use base 'Bio::Root::Root';

use File::Spec;
use File::Copy;

our $VERSION = '0.1';

sub is_integer {
    my ( $self, $variable, $desc) = @_;
    return $self->_check_variable( variable => $variable, desc => $desc, regex
    => 'int');
}

sub is_word {
    my ( $self, $variable, $desc) = @_;
    return $self->_check_variable( variable => $variable, desc => $desc, regex
    => 'word');
}

sub is_path {
    my ( $self, $variable, $desc) = @_;
    return $self->_check_variable( variable => $variable, desc => $desc, regex
    => 'path');
}

sub is_real {
    my ( $self, $variable, $desc) = @_;
    return $self->_check_variable( variable => $variable, desc => $desc, regex
    => 'real' );
}

sub is_sentence {
    my ( $self, $variable, $desc) = @_;
    return $self->_check_variable( variable => $variable, desc => $desc, regex
    => 'sentence' );
}

###############################################################################
# Usage      : $self->_check_variable(%args)
# Purpose    : make perl -T happy: get only "secure" symbols for system()
# Returns    : the specified variable if it is ok
# Parameters : ( variable => $v, regex => 'int', desc => $desc )
# Throws     : Bio::Root::BadParameter if variable is not ok
# Comments   : if variable is undef, then return value is undef
# See also   : http://gunther.web66.com/FAQS/taintmode.html

sub _check_variable {
   my ( $self, %args ) = @_;
   
   if ( !%args  || !defined $args{regex} ) {
      $self->throw( -class => 'Bio::Root::BadParameter',
                 -text  => 'Missing arguments: require hash with keys "regex",
                 "variable" and optional "desc"',);
   } 
   return if !defined $args{variable};

   if (!defined $args{desc}) {
      $args{desc} = 'Variable';
   }
   my $value;
   
   if ($args{regex} eq 'int') {
       ( $value ) = $args{variable} =~ m{ ( \A \d+ \z ) }xms;
   }    
   elsif ($args{regex} eq 'word') {
       ( $value ) = $args{variable} =~ m{ ( \A [\w.\-]+ \z ) }xms;
   }    
   elsif ($args{regex} eq 'path') {
       ( $value ) = $args{variable} =~ m{ ( \A [\w.\-/]+ \z ) }xms;
   }    
   elsif ($args{regex} eq 'sentence') {
       ( $value ) = $args{variable} =~ m{  \A ([\w.\-/|:(),;]+)  }xms;
   }    
   else {
       ( $value ) = $args{variable} =~ m{ ( \A -? \d+ \.? \d* \z ) }xms;
   }
   
#   my ( $value ) = $args{variable} =~ m{ (\d+) }xms;
   use Scalar::Util qw/tainted/;

   warn 'TAINTED' if tainted $args{regex};
   if ( !defined $value  || $value ne $args{variable}) {
      $value = '*undef*' if !defined $value;
      $self->throw( -class => 'Bio::Root::BadParameter',
                  -text  => $args{desc} . ' looks not valid.',
                  -value => $args{variable} . ' vs ' . $value );
   }              
   return $value;
}

1;# Magic true value required at end of module
__END__

=head1 NAME

Bio::Grep::Root - superclass for all Bio::Grep* packages


=head1 DESCRIPTION

This is the superclass adds some useful methods to packages.
Don't call them

=head1 METHODS

=over

=item C<is_integer($var_to_check, 'Description')>

Returns the variable, but it is now not tainted anymore. Throws an exception
if the specified variable is not an integer. If a second argument is passed,
it will be used in the exception text.

=item C<is_real($var_to_check, 'Description')>

Returns the variable, but it is now not tainted anymore. Throws an exception
if the specified variable is not a number. If a second argument is passed,
it will be used in the exception text.

=item C<is_word($var_to_check, 'Description')>

Returns the variable, but it is now not tainted anymore. Throws an exception
if the specified variable is not a word. If a second argument is passed,
it will be used in the exception text.

=item C<is_path($var_to_check, 'Description')>

Returns the variable, but it is now not tainted anymore. Throws an exception
if the specified variable is not a valid path (It is word but accepts also
'/'). If a second argument is passed, it will be used in the exception text.

=item C<is_sentence($var_to_check, 'Description')>

Returns the variable, but it is now not tainted anymore. Throws an exception
if the specified variable is not a valid sentence (It is word but accepts also
'/:,;|'). If a second argument is passed, it will be used in the exception text.

WARNING: make sure that command line arguments are quoted ( my $command = "...
 '$sentence' ") 

=back

=head1 INTERNAL METHODS

Only packages should call them

=over 

=item C<_check_variable( variable => $var_to_check, regex => $regex, desc =>
'Description')>

=back

=head1 DIAGNOSTICS

=over

=item C<Bio::Root::BadParameter>

Some variable does not match the specified regular expression.

=back

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
