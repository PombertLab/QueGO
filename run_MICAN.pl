#!/usr/bin/perl
## Pombert Lab 2022

my $name = "run_MICAN.pl";
my $version = "0.1.1";
my $updated = "2022-07-23";

use strict;
use warnings;
use PerlIO::gzip;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Calculates template model (TM) score with MICAN on structural
		matches identified with foldseek or GESAMT

USAGE

OPTIONS
-r (--results_dir)	RESULTS directory within STRUCTURAL_HOMOLOGY created by run_QueGO.pl
-u (--uniprot_pdb)	PDB directory withing UNIPROT_SCRAP_RESULTS created by run_QueGO.pl
-p (--predict_dir)	Directory(s) containing predicted structures

EXIT

my $results_dir;
my $uniprot_dir;
my @predicted_dirs;

GetOptions(
	'r|results_dir=s' => \$results_dir,
	'u|uniprot_pdb=s' => \$uniprot_dir,
	'p|predict_dir=s{1,}' => \@predicted_dirs,
);

my %predicted_dirs;
foreach my $dir (@predicted_dirs){
	my @data = split(/\//,$dir);
	$predicted_dirs{$data[-1]} = $dir;
}

unless (-d "temp"){
	make_path("temp",{mode=>0755});
}

opendir(ODIR,$results_dir) or die "Cannot access $results_dir: $!\n";
while (my $hom_tool_dir = readdir(ODIR)){
	if ((-d $results_dir."/".$hom_tool_dir) && ($hom_tool_dir !~ /^\./)){
		opendir(MDIR,$results_dir."/".$hom_tool_dir) or die "Cannot access $results_dir/$hom_tool_dir: $!\n";;
		while (my $structure_set_dir = readdir(MDIR)){
			if ((-d $results_dir."/".$hom_tool_dir."/".$structure_set_dir) && ($structure_set_dir !~ /^\./)){
				unless (-d $results_dir."/".$hom_tool_dir."_w_MICAN/".$structure_set_dir){
					make_path($results_dir."/".$hom_tool_dir."_w_MICAN/".$structure_set_dir,{mode=>0755});
				}
				opendir(IDIR,$results_dir."/".$hom_tool_dir."/".$structure_set_dir) or die "Cannot access $results_dir/$hom_tool_dir/$structure_set_dir";
				
				my $file_count = () = readdir(IDIR);
				rewinddir(IDIR);
				my $file_counter = 0;
				
				while (my $file = readdir(IDIR)){
					if ((-f $results_dir."/".$hom_tool_dir."/".$structure_set_dir."/".$file) && ($file ne "error.log")){
						my ($outfile) = $file =~ /^(\w+)/;
						$file_counter++;

						my $done = int(($file_counter/$file_count)*100);
						my $remaining = 100-$done;

						my $bar = "\t[".("|"x$done).("."x$remaining)."]\t($file_counter/$file_count)\n";

						print($bar."Working on $outfile...\n");

						unless (-f "$results_dir/${hom_tool_dir}_w_MICAN/$structure_set_dir/${outfile}_w_tmscore.fseek.gz"){
							my $gzip = "";
							if ($file =~ /\.gz$/){
								$gzip = ":gzip"
							}
							open IN, "<$gzip", $results_dir."/".$hom_tool_dir."/".$structure_set_dir."/".$file;
							my @results;
							while (my $line = <IN>){
								my $target_pdb = $uniprot_dir."/";
								my $pred_pdb = $predicted_dirs{$structure_set_dir}."/";
								chomp($line);
								if ($hom_tool_dir eq "FOLDSEEK"){
									my @data = split("\t",$line);

									$target_pdb .= ("/".$data[0]);
									my $temp_target = "temp/target.pdb";
									system "zcat -f $target_pdb > $temp_target";
									$pred_pdb .= ("/".$data[1]);
									my $temp_pred = "temp/pred.pdb";
									system "zcat -f $pred_pdb > $temp_pred";

									my $mican_result = `mican -s $temp_target $temp_pred -n 1`;
									
									my @mican_data = split("\n",$mican_result);

									my $grab;
									my $rank;
									my $sTMscore;
									my $TMscore;
									my $Dali_Z;
									my $SPscore;
									my $Length;
									my $RMSD;
									my $Seq_Id;
									foreach my $line (@mican_data){
										chomp($line);
										if ($line =~ /Rank\s+sTMscore/){
											$grab = 1;
										}
										if (($grab) && ($line =~ /^\s+(1.*)/)){
											undef($grab);
											($rank,$sTMscore,$TMscore,$Dali_Z,$SPscore,$Length,$RMSD,$Seq_Id) = split(/\s+/,$1);
											push(@data,$TMscore);
											push(@results,[@data]);
										}
									}
								}
							# 	elsif ($hom_tool_dir eq "GESAMT"){

							# 	}
							}
							close IN;
							open OUT, ">", "$results_dir/${hom_tool_dir}_w_MICAN/$structure_set_dir/${outfile}_w_tmscore.fseek";
							foreach my $line (sort{@{$b}[-1] <=> @{$a}[-1]}@results){
								print OUT (join("\t",@{$line})."\n");
							}
							close OUT;
							system "gzip $results_dir/${hom_tool_dir}_w_MICAN/$structure_set_dir/${outfile}_w_tmscore.fseek";
						}
					}
				}
				close IDIR;
			}
		}
		close MDIR;
	}
}
close ODIR;

system("rm -r temp");