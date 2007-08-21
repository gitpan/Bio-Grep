use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Grep',
    license             => 'perl',
    dist_author         => 'Markus Riester <mriester@gmx.de>',
    dist_version_from   => 'lib/Bio/Grep.pm',
    requires            => {
                            'Class::MethodMaker' => 2.02,
                            'version'            => 0,
    		                'IO::String'         => 1.0, 
                            'Getopt::Long'       => 2.32,   
                            'Pod::Usage'         => 0,   
                            'UNIVERSAL::require' => 0.11,   
                            'Carp::Assert'       => 0,
                            'ToolSet'            => 0.13,
                            'Bio::Root::Version' => 0,
                            },
    script_files        => [ 'script/bgrep' ], 
    add_to_cleanup      => [ 'Bio-Grep-*' ],
);

$builder->create_build_script();
