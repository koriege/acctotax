#! /usr/bin/perl

use strict;
use warnings;
use Tree::Simple;
use Bio::Tree::Node;
use Bio::TreeIO;

my $leafMap;
my $tree = Bio::Tree::Node->new();
open F , '</home/koriege/infatility_analysis/nr_blast_out.uniq.tab' or die $!;	
	while(<F>){
		my @line = split(/\s+/,$_);		
		next if $#line < 3;
		next if $line[3]=~/(C|c)hordata/;
		next if $line[3]=~/(A|a)rtificial/;
		next if $line[3]=~/(U|u)nclassified/;
		next if $line[2] > $ARGV[1];
		next if length $line[1] < $ARGV[0];
		
		my @l = split(/;/,$line[3]);
		for (1..$ARGV[2]){
			if($#l>2){
				pop @l;
			}
		}
		
		$_=~s/\W//g for @l;
		&insert3($tree, @l); 
	}			
close F;

my @removeNodes;
my $obj = Bio::Tree::Tree->new(-root => $tree);
for ($tree->get_all_Descendents){
	if ($_->is_Leaf){
		$_->id($_->id.'_#'.$leafMap->{$_->id});	
	}
}
$obj->contract_linear_paths();

#my $treeio = Bio::TreeIO->new(-format => 'svggraph' , -file => '>/home/koriege/infatility_analysis/'.$ARGV[0].'_'.$ARGV[1].'.svg');
#$treeio->write_tree($obj);
my $treeio = Bio::TreeIO->new(-format => 'newick' , -file => '>/home/koriege/infatility_analysis/trees/pop'.$ARGV[2].'_length'.$ARGV[0].'_cutoff'.$ARGV[1].'.nwk');
$treeio->write_tree($obj);

sub insert3 {
	my ($ref, $head, @tail) = @_;
	if ( @tail ) {
		my @tmp = grep {$_->id eq $head} $ref->each_Descendent;
		if($#tmp==-1){
			my $node = Bio::Tree::Node->new(-id => $head);
			$ref->add_Descendent($node);
			&insert3(  $node , @tail );			
		} else {
			&insert3(  $tmp[0] , @tail );
		} 
	}else{
		$leafMap->{$head}++;
		my @tmp = grep {$_->id eq $head} $ref->each_Descendent;
		$ref->add_Descendent(Bio::Tree::Node->new(-id => $head)) if $#tmp==-1;		
	}
}