use 5.008001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Bio::Grep',
    VERSION_FROM      => 'lib/Bio/Grep.pm', # finds $VERSION
    PREREQ_PM         => { 
                            'Class::MethodMaker' => 2.02,
                            'version'            => 0,
    		                'IO::String'         => 1.0, 
                            'Getopt::Long'       => 2.32,   
                            'Pod::Usage'         => 0,   
                          }, 
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Bio/Grep.pm', # retrieve abstract from module
       AUTHOR         => 'Markus Riester <mriester@gmx.de>') : ()),
    EXE_FILES         => [ 'script/bgrep' ], 
);
