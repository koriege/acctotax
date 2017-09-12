use ExtUtils::MakeMaker;
use strict;
use warnings;

my @scripts = qw (acctotax.pl);

WriteMakefile(
	'NAME' => 'AccToTax',
	'AUTHOR' => q{Konstantin Riege},
	'VERSION_FROM' => 'acctotax.pl',
	'LICENSE' => 'MIT',
	'PREREQ_PM' => {
		'Getopt::Long' => '0',
		'Bio::Root::Version' => '0',
		'Tree::Simple' => '0',
		'Bio::Tree::Node' => '0',
		'Bio::TreeIO' => '0',
		'Bio::DB::EUtilities' => '0',
		'Try::Tiny' => '0',
		'XML::Simple' => '0',
		'Getopt::Long' => '0',
		'Bio::DB::Taxonomy' => '0',
		'Bio::Tree::Tree' => '0',
		'File::Path' => '0',
		'File::Spec::Functions' => '0',
		'File::Basename' => '0',
	},
	'MIN_PERL_VERSION'	=>	'5.010',
	'EXE_FILES'	=>	[ map { "$_" } @scripts ],
	'CONFIGURE_REQUIRES'	=>	{
		'ExtUtils::MakeMaker'	=>	0
	},
	'BUILD_REQUIRES'	=>	{
		'Test::More' => 0
	},
	dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
	# clean => { FILES => '*' },
);