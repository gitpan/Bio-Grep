package Bio::Grep::Filter::FilterI;

use strict;
use warnings;

use version; our $VERSION = qv('0.5.0');

use base 'Bio::Root::Root';

use Class::MethodMaker
 [ new => 'new2',
   scalar => [qw / search_result message delete / ],
   hash => [qw / supports_alphabet/ ],
   abstract => [qw / filter new reset / ],
 ];


1;# Magic true value required at end of module
__END__


=head1 NAME

Bio::Grep::Filter::FilterI - Superclass for all filter modules   


=head1 DESCRIPTION

B<Bio::Grep::Filter::FilterI> is the superclass for all filter modules. 
Don't use this directly. 

Use this as interface for new filters (see L<Bio::Grep::Filter::MIRNAFilter> for
details). A Filter module implements a filter function that
returns 1 if query and subject pass the filter, 0 otherwise.
If the member variable "delete" is 1, then this subject won't
be in the search results of the back-end.


=head1 METHODS

=over

=item C<$filter-E<gt>delete()>

Get/set delete. If this is 0, then the search result will not be deleted  
(but you will still have the filter message). 
Otherwise, this hit won't be included in the search results of the back-end. 
Default ist 1.

=back


=head1 ABSTRACT METHODS

Every filter must implement this methods. 

=over

=item C<new()>

This function constructs a filter object.

=item C<$filter-E<gt>filter()>

This function returns 1 if query and subject pass the filter, 0 otherwise. You
have to set the search result with the function
C<$filter-E<gt>search_result> before. L<Bio::Grep::Backends::BackendI>
takes care of that.

Internal method: Only L<Bio::Grep::Backends::BackendI> should call this method.

=item C<$filter-E<gt>reset()>

Get/set reset. A flag needed by some Filters like FilterRemoveDuplicates to tell them, it
is a new search, forget everything. 

Internal method: Only L<Bio::Grep::Backends::BackendI> should call this method.

=item C<$filter-E<gt>supports_alphabet()>

Get supported alphabets. Returns a hash. Keys are the supported alphabets.

    my $can_filter_proteins = $filter->supports_alphabet_exists('protein');

=back

=head1 INTERNAL METHODS

Only L<Bio::Grep::Backends::BackendI> should call them directly.

=over

=item C<$filter-E<gt>message()>

Get/set the message. This is a string with the reason for rejecting the search
result. 

=item C<$filter-E<gt>search_result()>

Get/set the search_result. A L<Bio::Grep::Container::SearchResult>
object. 

=back

=head1 SEE ALSO

L<Bio::Grep::Backends::BackendI>
L<Bio::Grep::Container::SearchResult>

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
