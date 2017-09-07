#! /usr/bin/env perl
use v5.10;
use strict;
use warnings;
use Tree::Simple;
use Bio::Tree::Node;
use Bio::TreeIO;
use Bio::DB::EUtilities;
use Try::Tiny;
use XML::Simple;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;

my ($protids,$names,$taxids,$nucleotides,$gids,$out,$taxdb,$nodesfile,$namesfile,$help) = ('','','','','','');
(Getopt::Long::Parser->new)->getoptions (
		'dbnodes|ncbinodes=s' => \$nodesfile,
		'dbnames|ncbinames=s' => \$namesfile,
		'n|names=s' => \$names,
		't|taxids=s' => \$taxids,
		'n|nucleotides=s' => \$nucleotides,
		'g|genomeids=s' => \$gids,
		'p|proteins=s' => \$protids,
		'o|out=s' => \$out,
		'h|help' => \$help
) or do {print "Usage: -out /my/tree/file -[names|taxids|nucleotides|genomeids|proteins] <\'list,of,values\'>\noptional - use local taxonomy files: -dbnodes nodes.dmp -dbnames names.dmp\n"; exit};
do {print "Usage: -out /my/tree/file -[names|taxids|nucleotides|genomeids|proteins] <\'list,of,values\'>\noptional - use local taxonomy files: -dbnodes nodes.dmp -dbnames names.dmp\n"; exit} if $help || ! ($names || $taxids || $nucleotides || $gids || $protids);

if ($nodesfile && $namesfile){
	mkdir "/tmp/$$";
	$taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namesfile, -directory => "/tmp/$$" , -force => 1 , -verbose => -1 ) ;
}

my $leafMap;
my $tree = Bio::Tree::Node->new();

my $i=0;
my @queries = split(/\s*,\s*/,$gids);
try {
	for (@queries){	
		#print ''.(++$i)." of $#queries\n";
		print $_.":".&getNameFromID($_)." " for &getLineageNodes( &getIDfromGenomeID( $_ ) );
		say "";
		# &insert($tree,&getLineageNodes( &getIDfromGenomeID( $_ ) ) );
	}	
	$i=0;
	@queries = split(/\s*,\s*/,$nucleotides);
	for (@queries){	
		#print ''.(++$i)." of $#queries\n";
		print $_.":".&getNameFromID($_)." " for &getLineageNodes( &getIDfromNucleotide( $_ ) );
		say "";
		# &insert($tree,&getLineageNodes( &getIDfromAccession( $_ ) ) );
	}
	$i=0;
	@queries = split(/\s*,\s*/,$names);
	for (@queries){	
		#print ''.(++$i)." of $#queries\n";
		print $_.":".&getNameFromID($_)." " for &getLineageNodes( &getIDfromName( $_ ) );
		say "";
		# &insert($tree,&getLineageNodes( &getIDfromName( $_ ) ) );
	}
	$i=0;
	@queries = split(/\s*,\s*/,$taxids);
	for (@queries){	
		# print ''.(++$i)." of $#queries\n";
		print $_.":".&getNameFromID($_)." " for &getLineageNodes( $_ );
		say "";
		# &insert($tree,&getLineageNodes( $_ ) );
	}
	$i=0;
	@queries = split(/\s*,\s*/,$protids);
	for (@queries){	
		# print ''.(++$i)." of $#queries\n";
		print $_.":".&getNameFromID($_)." " for &getLineageNodes( &getIDfromProtein( $_ ) );
		say "";
		# &insert($tree,&getLineageNodes( $_ ) );
	}
} catch {
	say "";
};

# my $obj = Bio::Tree::Tree->new(-root => $tree);
# $obj->contract_linear_paths();
# for ($tree->get_all_Descendents){
# 	#if ($_->is_Leaf){
# 		$_->id($_->id.'_'.&getNameFromID($_->id));
# 	#}
# }

# my $treeio = Bio::TreeIO->new(-format => 'newick' , -file => '>'.$out);
# $treeio->write_tree($obj);

unlink "/tmp/$$";

sub insert {
	my ($ref, $head, @tail) = @_;

	return unless $head;

	if ( @tail ) {
		my @tmp = grep {$_->id eq $head} $ref->each_Descendent;
		if($#tmp==-1){
			my $node = Bio::Tree::Node->new(-id => $head);
			$ref->add_Descendent($node);
			&insert(  $node , @tail );			
		} else {
			&insert(  $tmp[0] , @tail );
		} 
	}else{
		$leafMap->{$head}++;
		my @tmp = grep {$_->id eq $head} $ref->each_Descendent;
		$ref->add_Descendent(Bio::Tree::Node->new(-id => $head)) if $#tmp==-1;		
	}
}

sub getLineageNodes {
	my ($taxid) = @_;

	return () unless $taxid;

	my $taxon;
	$taxon = $taxdb->get_taxon(-taxonid => $taxid) if $taxdb;
	my @nodes;
	if ($taxon){
		my $tree_functions = Bio::Tree::Tree->new( -verbose => -1);
		@nodes = $tree_functions->get_lineage_nodes($taxon);
		push @nodes, $taxid;
	}
	if ($#nodes == -1){
		my $errorCounter=0;
		my $error=1;
		while($error && $errorCounter < 10){
			$errorCounter++;
			$error=0;
			try{			
				my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch', -email => 'mymail@foo.bar', -db => 'taxonomy', -id => $taxid  , -verbose => -1);
				my $res = $factory->get_Response->content;
				my $data = XMLin($res);
				do {push @nodes, $_->{TaxId} for @{$data->{Taxon}->{LineageEx}->{Taxon}} } if ref $data;
			} catch {
				$error=1;
			};
		}
		push @nodes, $taxid;
	}
	
	return @nodes;
}

sub getIDfromNucleotide {
	my ($acc) = @_;	

	my $error=1;
	my $errorCounter=0;
	my $taxid;
	while($error && $errorCounter < 10){
		$errorCounter++;
		$error=0;
		try{
			my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch', -email => 'mymail@foo.bar', -db => 'nuccore', -term => $acc  , -verbose => -1);
			my @uids = $factory->get_ids;
			$factory->reset_parameters(-eutil => 'esummary', -email => 'mymail@foo.bar', -db => 'nuccore', -id => \@uids);	
			my $docsum = $factory->next_DocSum;
			($taxid) = $docsum->get_contents_by_name('TaxId');								

		} catch {					
			$error=1;
		};
	}
	
	return $taxid;
}

sub getIDfromProtein {
	my ($acc) = @_;	

	my $error=1;
	my $errorCounter=0;
	my $taxid;
	while($error && $errorCounter < 10){
		$errorCounter++;
		$error=0;
		try{
			my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch', -email => 'mymail@foo.bar', -db => 'protein', -term => $acc  , -verbose => -1);
			my @uids = $factory->get_ids;
			$factory->reset_parameters(-eutil => 'esummary', -email => 'mymail@foo.bar', -db => 'protein', -id => \@uids);	
			my $docsum = $factory->next_DocSum;
			($taxid) = $docsum->get_contents_by_name('TaxId');
		} catch {					
			$error=1;
		};
	}
	
	return $taxid;
}

sub getNameFromID {
	my ($taxid) = @_;	

	my $taxon;
	$taxon = $taxdb->get_taxon(-taxonid => $taxid) if $taxdb;
	my $name;
	if ($taxon) {
		$name = $taxon->scientific_name;
	} else {		
		my $errorCounter=0;
		my $error=1;		
		while($error && $errorCounter < 10){
			$errorCounter++;
			$error=0;
			try{
				my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary', -email => 'mymail@foo.bar', -db => 'taxonomy', -id => $taxid  , -verbose => -1);	
				($name) = $factory->next_DocSum->get_contents_by_name('ScientificName');
			} catch {					
				$error=1;
			};
		} 		
	}
	if ($name){
		$name=~s/\s+/_/g;
		$name=~s/\.\./\./g;
		$name=~ s/[^a-zA-Z0-9_]*//g;		
		$name=ucfirst($name);
		return $name;
	} else {
		return $taxid;
	}		
					
}

sub getIDfromGenomeID {
	my ($acc) = @_;	

	my $error=1;
	my $errorCounter=0;
	my $taxid;
	while($error && $errorCounter < 10){
		$errorCounter++;
		$error=0;
		try{
			my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary', -email => 'mymail@foo.bar', -db => 'genome', -id => $acc  , -verbose => -1);
			# my @uids = $factory->get_ids;
			# $factory->reset_parameters(-eutil => 'esummary', -email => 'mymail@foo.bar', -db => 'nuccore', -id => \@uids);	
			my $docsum = $factory->next_DocSum;
			# print $docsum;
			($taxid) = &getIDfromAccession($docsum->get_contents_by_name('Assembly_Accession'));	
		} catch {					
			$error=1;
		};
	}	
	
	return $taxid;
}

sub getIDfromName {
	my ($query) = @_;

	my $taxid;
	my $errorCounter=0;
	my $error=1;
	while($error && $errorCounter < 10){
		$errorCounter++;
		$error=0;
		try{
			my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch', -db => 'taxonomy', -email => 'mymail@foo.bar', -term  => $query , -verbose => -1);
			#todo
			#handle selection if mistypos								
			for($factory->get_ids){						
				my @nodes = &getLineageNodes($_);					
				next if $#nodes<1;
				$taxid = $_;
				last if $taxid;
			}												
		} catch {					
			$error=1;
		};
	} 

	return $taxid;
}