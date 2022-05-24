#!/usr/bin/perl
## Pombert Lab 2022

my $name = "parse_GESAMT_results.pl";
my $version = "0.1";
my $updated = "2022-02-20";

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	

USAGE		${name} \\
			  -r Queri3D/RESULTS/ALPHAFOLD Queri3D/RESULTS/RAPTORX

OPTIONS
-r (--results)		Directory(s) containing GESAMT results
-q (--qscore)		Q-score cut-off value [Default = 0.3]
-b (--best)		Keep best 'X' results [Default = 25]
-a (--all)		Keep all results
-o (--outdir)		Output directory [Default = "./GESAMT_Results_Parsed"]
EXIT

die("\n$usage\n") unless(@ARGV);

my @results;
my $qscore_cut = 0.3;
my $best = 25;
my $all;
my $outdir = "GESAMT_Results_Parsed";

GetOptions(
	'r|results=s{1,}' => \@results,
	'q|qscore=s' => \$qscore_cut,
	'b|best=s' => \$best,
	'a|all' => \$all,
	'o|outdir=s' => \$outdir,
);

unless(-d $outdir){
	make_path($outdir,{mode=>0755}) or die("Unable to create directory $outdir: $!\n");
}

my %results;
foreach my $directory (@results){
	opendir(DIR,$directory);
	my @source = split(/\//,$directory);
	my $source = $source[-1];
	while(my $file = readdir(DIR)){
		unless(-d $file){
			my ($query) = $file =~ /(\w+).normal.gesamt/;
			open IN, "<", "$directory/$file";
			while(my $line = <IN>){
				chomp($line);
				unless(($line =~ /^\#/)||($line eq '')){
					# print($line."\n");
					my @data = split('\s+',$line);
					my $qscore;
					my $rmsd;
					my $seq_id;
					my $n_align;
					my $nRes;
					my $source_file;
					if(scalar(@data) == 8){
						(undef,undef,$qscore,$rmsd,$seq_id,$n_align,$nRes,$source_file) = @data;
					}
					if(scalar(@data) == 9){
						(undef,undef,undef,$qscore,$rmsd,$seq_id,$n_align,$nRes,$source_file) = @data;
					}
					if(scalar(@data) == 10){
						(undef,undef,undef,undef,$qscore,$rmsd,$seq_id,$n_align,$nRes,$source_file) = @data;
					}
					my ($predicted) = $source_file =~ /^(\w+)/;
					if($qscore >= $qscore_cut){
						push(@{$results{$query}{$predicted}},($source,$qscore,$rmsd,$seq_id,$n_align,$nRes));
					}
				}
			}
		}
	}
}

open ALL, ">", "$outdir/All_Parsed_Results.matches";
print ALL ("### Locus\tSource\tQ-Score\tr.m.s.d\tSeq. Id.\tNalign\tnRes\n\n");
foreach my $query (sort(keys(%results))){
	print ALL ("## $query\n");
	my $query_count = 0;
	foreach my $predicted (sort{$results{$query}{$b}[1] <=> $results{$query}{$a}[1]}(keys(%{$results{$query}}))){
		my ($source,$qscore,$rmsd,$seq_id,$n_align,$nRes) = @{$results{$query}{$predicted}};
		if($all){
			print ALL ("$predicted\t$source\t$qscore\t$rmsd\t$seq_id\t$n_align\t$nRes\n");
		}
		elsif($query_count < $best){
			print ALL ("$predicted\t$source\t$qscore\t$rmsd\t$seq_id\t$n_align\t$nRes\n");
		}
		$query_count++;
	}
	print ALL ("\n");
}
close ALL;