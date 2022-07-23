#!/usr/bin/perl
## Pombert Lab 2022

my $name = "parse_3D_homology_results.pl";
my $version = "0.2.1";
my $updated = "2022-07-22";

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use PerlIO::gzip;

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	

USAGE		${name} \\
			  -g Queri3D/RESULTS/ALPHAFOLD Queri3D/RESULTS/RAPTORX

OPTIONS
-g (--gesamt)		Directory(s) containing GESAMT results
-f (--foldseek)		Directory(s) containing FoldSeek results
-q (--qscore)		Q-score cut-off value [Default = 0.3]
-t (--tm)		TM cut-off value [Default = 0.3]
-b (--best)		Keep best 'X' results [Default = 25]
-a (--all)		Keep all results
-o (--outdir)		Output directory [Default = "./3D_homology_results_parsed"]
EXIT

die("\n$usage\n") unless(@ARGV);

my %result_dirs;
my $qscore_cut = 0.3;
my $tm_cut = 0.3;
my $best = 25;
my $all;
my $outdir = "3D_homology_results_parsed";

GetOptions(
	'g|gesamt=s{1,}' => \@{$result_dirs{"GESAMT"}},
	'f|foldseek=s{1,}' => \@{$result_dirs{"FoldSeek"}},
	'q|qscore=s' => \$qscore_cut,
	'e|eval=s' => \$tm_cut,
	'b|best=s' => \$best,
	'a|all' => \$all,
	'o|outdir=s' => \$outdir,
);

unless(-d $outdir){
	make_path($outdir,{mode=>0755}) or die("Unable to create directory $outdir: $!\n");
}

my %results;
foreach my $predictor (keys(%result_dirs)){
	if ($predictor eq "GESAMT"){
		foreach my $directory (@{$result_dirs{$predictor}}){
			if (-d $directory){
				opendir(DIR,$directory);
				my @source = split(/\//,$directory);
				my $pred_struct_source = $source[-1];

				while(my $file = readdir(DIR)){
					
					unless(-d $file){

						my $query_struct;
						my $gzip;

						if ($file =~ /(\w+).normal.gesamt.gz$/){
							$query_struct = $1;
							$gzip = ":gzip";
						}
						elsif ($file =~ /(\w+).normal.gesamt$/){
							$query_struct = $1;
						}

						open IN, "<", "$directory/$file";
						while(my $line = <IN>){
							chomp($line);
							unless(($line =~ /^\#/)||($line eq '')){
								# print($line."\n");
								my @data = split('\s+',$line);

								my ($qscore,$rmsd,$seq_id,$n_align,$nRes,$predicted_file) = @data[($#data-5)..$#data];

								my ($predicted_structure) = $predicted_file =~ /^(\w+)/;
								if($qscore >= $qscore_cut){
									push(@{$results{$predictor}{$query_struct}{$predicted_structure}},($pred_struct_source,$qscore,$rmsd,$seq_id,$n_align,$nRes));
								}
							}
						}
					}
				}
			}
		}
	}
	elsif($predictor eq "FoldSeek"){
		foreach my $directory (@{$result_dirs{$predictor}}){
			if (-d $directory){
				opendir(DIR,$directory);
				my @source = split(/\//,$directory);
				my $pred_struct_source = $source[-1];
				while(my $file = readdir(DIR)){
					unless(-d $file || ($file eq "error.log")){

						my $query_struct;
						my $gzip;

						if ($file =~ /(\w+)_w_MICAN.fseek.gz$/){
							$query_struct = $1;
							$gzip = ":gzip";
						}
						elsif ($file =~ /(\w+)_w_MICAN.fseek$/){
							$query_struct = $1;
						}

						open IN, "<$gzip", "$directory/$file";
						while(my $line = <IN>){
							chomp($line);
							unless(($line =~ /^\#/)||($line eq '')){
								my @data = split('\s+',$line);
								my ($predicted_structure) = $data[1] =~ /^(\w+)/;
								if ($data[12] <= $tm_cut){
									push(@{$results{$predictor}{$query_struct}{$predicted_structure}},($pred_struct_source,@data[2..12]));
								}
							}
						}
					}
				}
			}
		}
	}
}

foreach my $predictor (keys(%results)){
	open ALL, ">", "$outdir/${predictor}_parsed_results.matches";
	if ($predictor eq "GESAMT"){
		print ALL ("### Locus\tSource\tQ-Score\tr.m.s.d\tSeq. Id.\tNalign\tnRes\n\n");
		foreach my $query_struct (sort(keys(%{$results{$predictor}}))){
			print ALL ("## $query_struct\n");
			my $query_count = 0;
			foreach my $predicted_structure (sort{$results{$predictor}{$query_struct}{$b}[1] <=> $results{$predictor}{$query_struct}{$a}[1]}(keys(%{$results{$predictor}{$query_struct}}))){
				my ($pred_struct_source,$qscore,$rmsd,$seq_id,$n_align,$nRes) = @{$results{$predictor}{$query_struct}{$predicted_structure}};
				if($all){
					print ALL ("$predicted_structure\t$pred_struct_source\t$qscore\t$rmsd\t$seq_id\t$n_align\t$nRes\n");
				}
				elsif($query_count < $best){
					print ALL ("$predicted_structure\t$pred_struct_source\t$qscore\t$rmsd\t$seq_id\t$n_align\t$nRes\n");
				}
				$query_count++;
			}
			print ALL ("\n");
		}
	}
	elsif($predictor eq "FoldSeek"){
		print ALL ("### Locus\tSource\tfident\talnlen\tmismatch.\tgapopen\tqstart\tqend\ttstart\ttend\teval\tbits\ttmscore\n\n");
		foreach my $query_struct (sort(keys(%{$results{$predictor}}))){
			print ALL ("## $query_struct\n");
			my $query_count = 0;
			foreach my $predicted_structure (sort{$results{$predictor}{$query_struct}{$b}[11] <=> $results{$predictor}{$query_struct}{$a}[11]}(keys(%{$results{$predictor}{$query_struct}}))){
				my ($pred_struct_source,$fident,$alnlen,$mismatch,$gapopen,$qstart,$qend,$tstart,$tend,$eval,$bits,$tmscore) = @{$results{$predictor}{$query_struct}{$predicted_structure}};
				if($all){
					print ALL ("$predicted_structure\t$pred_struct_source\t$fident\t$alnlen\t$mismatch\t$gapopen\t$qstart\t$qend\t$tstart\t$tend\t$eval\t$bits\t$tmscore\n");
				}
				elsif($query_count < $best){
					print ALL ("$predicted_structure\t$pred_struct_source\t$fident\t$alnlen\t$mismatch\t$gapopen\t$qstart\t$qend\t$tstart\t$tend\t$eval\t$bits\t$tmscore\n");
				}
				$query_count++;
			}
			print ALL ("\n");
		}
	}
	close ALL;
}