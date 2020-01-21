#! /usr/bin/env perl
# (c) Konstantin Riege
use v5.10;
use strict;
use warnings;
use version; our $VERSION = qv('0.1.2');
use sigtrap qw(handler cleanup INT TERM KILL);
use Tree::Simple;
use Bio::Tree::Node;
use Bio::TreeIO;
use Bio::DB::EUtilities;
use Try::Tiny;
use XML::Simple;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);
use File::Basename;

sub usage {
	my $s = basename($0);
print <<EOF;
DESCRIPTION
AccToTax is able to infer lineages, a taxonomic tree and scientific names from the NCBI database using its web API, 
given a list of fuzzy names, nucleotide accession numbers, protein accession numbers and/or taxonomic IDs.

VERSION
$VERSION

SYNOPSIS 
$s -nam [name] -nam [name] -nuc [accNr] -tree [path] [parameter] > [results]

GENERAL OPTIONS
-h       | --help                      : this message
-tmp     | --tmpdir [path]             : optional - set path to a directory - default: ./
-dbnodes | --dbncbinodes               : optional - set path to nodesfile of a local NCBI taxonomy for speedup
-dbnames | --dbncbinames               : optional - set path to namesfile of a local NCBI taxonomy for speedup
-r       | --reconnects                : optional - set number of tries to connect to NCBI web API - default: 2

INPUT OPTIONS
-nam     | --name [name]               : (multiple) input of single, scientific (fuzzy) organism names
-namf    | --namesfile [path]          : input a file containing names (one per line)
-nuc     | --nucleotidesacc [accNr]    : (multiple) input of single nucleotide accession numbers
-nucf    | --nucleotideaccsfile [path] : input a file containing nucleotide accession numbers (one per line)
-pro     | --proteinacc [accNr]        : (multiple) input of single protein accession numbers
-prof    | --proteinaccsfile [path]    : input a file containing protein accession numbers (one per line)
-tax     | --taxid [taxid]             : (multiple) input of single taxonomic IDs
-taxf    | --taxidsfile [path]         : input a file containing taxonomic IDs (one per line)

OUTPUT OPTIONS
-on      | --outnames                  : optional - output also scientific names
-tree    | --treefile [path]           : optional - create a subtree of NCBI taxonomy and output newick format
-c       | --contract                  : optional - contract linear paths in subtree output

REFERENCES
(c) Konstantin Riege
konstantin{a}bioinf{.}uni-leipzig{.}de
EOF
	exit 0;
}

usage if $#ARGV == -1;
my ($tries, $tmp, $tree) = (2, 'tmp'.$$.int(rand(999999)), Bio::Tree::Node->new());
my ($tmpdir, $out, $taxdb, $contract, $onames, $nodesfile, $namesfile, 
	%idToName, %idToLineage,
	@names, @taxids, @nuclaccs, @protaccs);

(Getopt::Long::Parser->new)->getoptions(
	'h|help' => sub{&usage},
	'r|reconnects=i' => \$tries,
	'tmp|tmpdir=s' => \$tmpdir,
	'c|contract' => \$contract,
	'tree|treefile=s' => \$out,
	'on|outnames' => \$onames,
	'dbnodes|dbncbinodes=s' => \$nodesfile,
	'dbnames|dbncbinames=s' => \$namesfile,
	'nam|name=s' => \@names,
	'namf|namesfile' => sub {push @names, $_ while <>},
	'tax|taxid=s' => \@taxids,
	'taxf|taxidsfile=s' => sub {push @taxids, $_ while <>},
	'nuc|nucleotideacc=s' => \@nuclaccs,
	'nucf|nucleotideaccsfile=s' => sub {push @nuclaccs, $_ while <>},
	'pro|proteinacc=s' => \@protaccs,
	'prof|proteinaccsfile=s' => sub {push @protaccs, $_ while <>},
) or &cleanup;

try {
	if ($nodesfile && $namesfile){
		try { # try to set up local db objects
			$tmp = catdir($tmpdir,$tmp) if $tmpdir;
			make_path($tmp);
			$taxdb = Bio::DB::Taxonomy->new(
				-source => "flatfile", 
				-nodesfile => $nodesfile, 
				-namesfile => $namesfile, 
				-directory => $tmp, 
				-force => 1, 
				-verbose => -1
			);
		} catch {
			print STDERR ":WARNING: Failed to load local taxonomy database\n";
		};
	} else {
		print STDERR ":INFO: No local taxonomy database available\n";
	}

	# print lineage from taxid - i.e. use code reference to obtain taxid from input
	# thereby catch errors and print warnings
	&queryout(\&getIDfromName, $_) for @names;
	&queryout(\&getIDfromNuclacc, $_) for @nuclaccs;
	&queryout(\&getIDfromProtacc, $_) for @protaccs;
	&queryout(sub {return $_[0]}, $_) for @taxids;

	if ($out){ # prepare tree output
		make_path(dirname($tmp));
		my $obj = Bio::Tree::Tree->new(-root => $tree);
		$obj->contract_linear_paths() if $contract; # contract tree
		if ($onames){ # add names to ids
			for ($tree->get_all_Descendents){
				$_->id($_->id.':'.$idToName{$_->id});
			}
		}
		(Bio::TreeIO->new(-format => "newick" , -file => ">$out"))->write_tree($obj);
	}
	remove_tree($tmp);
} catch {
	print STDERR ":ERROR: $_";
	&cleanup;
};

# trap signals
sub cleanup {
	print STDERR ":ERROR: Abort\n";
	remove_tree($tmp);
	exit 1;
}

sub queryout {
	my ($fun, $query) = @_;
	try {
		my @nodes = &getLineageIDs(&$fun($query));
		for (@nodes){
			print $_;
			if ($onames){ # add names to output
				my $n;
				try {
					$n = &getNameFromID($_);
				} catch {
					print STDERR ":WARNING: $_";
					$n = "?";
				};
				$idToName{$_} = $n;
				print ":$n";
			}
			print " ";
		}
		print "\n";
		&treeinsert($tree, @nodes) if $out;
	} catch {
		print STDERR ":WARNING: $_";
	};
}

# recursively add a linage into Bio::Tree object
sub treeinsert {
	my ($ref, $head, @tail) = @_;
	return unless $head;

	if (@tail){
		my @tmp = grep {$_->id eq $head} $ref->each_Descendent;
		if($#tmp == -1){
			my $node = Bio::Tree::Node->new(-id => $head);
			$ref->add_Descendent($node);
			&treeinsert($node, @tail);
		} else {
			&treeinsert($tmp[0], @tail);
		} 
	} else {
		my @tmp = grep {$_->id eq $head} $ref->each_Descendent;
		$ref->add_Descendent(Bio::Tree::Node->new(-id => $head)) if $#tmp == -1;
	}
}

# remove special characters from scientific name output
sub flatname {
	my ($name) = @_;

	$name=~s/(^\s+|\s+$)//g;
	$name=~s/\s+/_/g;
	$name=~s/\.\./\./g;
	$name=~s/__/_/g;
	$name=~s/[^a-zA-Z0-9_]*//g;
	$name=ucfirst($name);
	return $name;
}

sub getNameFromID {
	my ($taxid) = @_;	
	
	my $name = $idToName{$taxid};
	return $name if $name;
	if ($taxdb) { # check this first
		my $taxon = $taxdb->get_taxon(-taxonid => $taxid);
		$name = $taxon->scientific_name if $taxon;
	}
	unless ($name){
		my ($try, $error) = (0, 1);
		while($error && (++$try) < $tries){
			$error = 0;
			try { # try to reach web Aapi
				my $factory = Bio::DB::EUtilities->new(
					-eutil => 'esummary', 
					-email => 'mymail@foo.bar', 
					-db => 'taxonomy', 
					-id => $taxid, 
					-verbose => -1
				);
				($name) = $factory->next_DocSum->get_contents_by_name('ScientificName');
			} catch {
                print STDERR ":ERROR: $_";
				$error = 1;
			};
		}
	}
	if ($name){
		return &flatname($name);
	} else {
		die "No name found for $taxid"; # above: catch error and print as warning
	}
}

sub getIDfromName {
	my ($query) = @_;

	my $taxid;
	if ($taxdb) { # check this first
		my $taxon = $taxdb->get_taxon(-name => $query);
		$taxid = $taxon->id if $taxon;
	}
	unless ($taxid){
		my ($try, $error) = (0, 1);
		while ($error && (++$try) < $tries){
			$error = 0;
			try { # try to reach web api
				my $factory = Bio::DB::EUtilities->new(
					-eutil => 'esearch', 
					-db => 'taxonomy', 
					-email => 'mymail@foo.bar', 
					-term => $query, 
					-verbose => -1
				);
				for ($factory->get_ids){
					my @nodes = &getLineageIDs($_);
					next if $#nodes < 1;
					$taxid = $_;
					push @{$idToLineage{$taxid}}, @nodes;
					last;
				}
			} catch {
                print STDERR ":ERROR: $_";
				$error = 1;
			};
		}
	}
	if ($taxid){
		return $taxid;
	} else {
		die "No taxid found for $query"; # above: catch error and print as warning
	}
}

sub getIDfromNuclacc { # web api only
	my ($query) = @_;

	my ($taxid, $name);
	my ($try, $error) = (0, 1);
	while ($error && (++$try) < $tries){
		$error = 0;
		try {
			my $factory = Bio::DB::EUtilities->new(
				-eutil => 'esearch', 
				-email => 'mymail@foo.bar', 
				-db => 'nuccore', 
				-term => $query, 
				-verbose => -1
			);
			my @uids = $factory->get_ids;
			$factory->reset_parameters(
				-eutil => 'esummary', 
				-email => 'mymail@foo.bar', 
				-db => 'nuccore', 
				-id => \@uids
			);
			while (my $docsum = $factory->next_DocSum) {
				my ($id) = $docsum->get_contents_by_name('TaxId');
				# check for valid lineage of taxid from query matching api returns
				my @nodes = &getLineageIDs($id); 
				next if $#nodes < 1;
				$taxid = $id;
				# store valid lineage for later reuse
				push @{$idToLineage{$taxid}}, @nodes;
				last;
				# $name = &flatname($docsum->get_contents_by_name('Title'));
			}
		} catch {
            print STDERR ":ERROR: $_";
			$error = 1;
		};
	}
	if ($taxid){
		# $idToName{$taxid} = $name if $name;
		return $taxid;
	} else {
		die "No taxid found for $query"; # above: catch error and print as warning
	}
}

sub getIDfromProtacc { # web api only
	my ($query) = @_;

	my ($taxid, $name);
	my ($try, $error) = (0, 1);
	while ($error && (++$try) < $tries){
		$error = 0;
		try {
			my $factory = Bio::DB::EUtilities->new(
				-eutil => 'esearch', 
				-email => 'mymail@foo.bar', 
				-db => 'protein', 
				-term => $query, 
				-verbose => -1
			);
			my @uids = $factory->get_ids;
			$factory->reset_parameters(-eutil => 'esummary', -email => 'mymail@foo.bar', -db => 'protein', -id => \@uids);	
			while (my $docsum = $factory->next_DocSum) {
				my ($id) = $docsum->get_contents_by_name('TaxId');
				# check for valid lineage of taxid from query matching api returns
				my @nodes = &getLineageIDs($id);
				next if $#nodes < 1;
				$taxid = $id;
				# store valid lineage for later reuse
				push @{$idToLineage{$taxid}}, @nodes;
				last;
				# $name = &flatname($docsum->get_contents_by_name('Title'));
			}
		} catch {
            print STDERR ":ERROR: $_";
			$error=1;
		};
	}
	if ($taxid){
		# $idToName{$taxid} = $name if $name;
		return $taxid;
	} else {
		die "No taxid found for $query"; # above: catch error and print as warning
	}
}

sub getLineageIDs {
	my ($taxid) = @_;

	my ($taxon, @nodes);
	if (exists $idToLineage{$taxid}){ # check for reuse of prior crated valid lineages
		return @{$idToLineage{$taxid}};
	}

	$taxon = $taxdb->get_taxon(-taxonid => $taxid) if $taxdb;
	if ($taxon){ # check this first
		my $tree_functions = Bio::Tree::Tree->new(-verbose => -1);
		@nodes = $tree_functions->get_lineage_nodes($taxon);
		push @nodes, $taxid;
	}
	if ($#nodes == -1){ # try to reach web api
		my ($error, $try) = (1, 0);
		while($error && (++$try) < $tries){
			$error = 0;
			try{
				my $factory = Bio::DB::EUtilities->new(
					-eutil => 'efetch', 
					-email => 'mymail@foo.bar', 
					-db => 'taxonomy', 
					-id => $taxid,
					-verbose => -1
				);
				my $res = $factory->get_Response->content;
				# parse xml return into hash
				my $data = XMLin($res);
				if (ref $data){
					for (@{$data->{Taxon}->{LineageEx}->{Taxon}}){
						push @nodes, $_->{TaxId};
						$idToName{$_->{TaxId}} = &flatname($_->{ScientificName});
					}
				}
			} catch {
                print STDERR ":ERROR: $_";
				$error = 1;
			};
		}
		push @nodes, $taxid;
	}
	if ($#nodes > -1){
		push @{$idToLineage{$taxid}}, @nodes;
		return @nodes;
	} else {
		die "No lineage found for $taxid"; # above: catch error and print as warning
	}
}

exit 0
