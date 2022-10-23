#!/usr/bin/perl

my $name = "create_visualizations.pl";
my $version = "0.7.6";
my $updated = "2022-09-03";

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd qw(abs_path);
use File::Path qw(make_path);

my $usage = << "EXIT";
NAME	${name}
VERSION	${version}
UPDATED	${updated}

SYNOPSIS	Creates ChimeraX visual comparisons between .pdb structures
		and matches found with QueGO

USAGE	${name} \\
		-m RESULTS/compiled_results.tsv \\
		-u /media/FatCat_1/julian/QueGO/UniProt \\
		-p /media/FatCat_1/Microsporidia/Intestinalis/3DFI/FOLDING/ALPHAFOLD_3D_PARSED

OPTIONS
-m (--match)	Foldseek/GESAMT match file parsed by organize_results.pl
-p (--prov)	Directory containing provided .pdb files
-u (--uni)	Directory containing UniProt scrap .pdb files
-o (--outdir)	Output directory for ChimeraX sessions [Default: ./3D_Visualizations]
EXIT
die "\n\n$usage\n\n" unless @ARGV;

my $match_file;
my @provided_struct;
my $uniprot_struct;
my $outdir = './3D_Visualizations';

GetOptions(
	'm|match=s' => \$match_file,
	'p|prov=s@{1,}' => \@provided_struct,
	'u|uni=s' => \$uniprot_struct,
	'o|out=s' => \$outdir,
);

my ($filename,$dir) = fileparse($0);
my $script = "$dir/chimerax_session_creator.py";

unless (-d $outdir){
	make_path($outdir,{mode => 0755}) or die "\n[ERROR]\tUnable to create $outdir: $!\n";
}

my %provided_struct;
foreach my $provided (@provided_struct){
	my ($dirname) = $provided =~ /\/(\w+)\/?$/;
	$provided_struct{$dirname} = $provided;
}

##
my %results;
open IN, "<", $match_file or die "Unable to read from $match_file: $!\n";
my $prot_name;
while (my $line = <IN>){

	chomp($line);
	if ($line =~ /^## (.*)\t/){
		$prot_name = $1;
		$prot_name =~ tr/ .,'"()/_/;
		make_path("$outdir/$prot_name",{mode => 0755});
	}
	elsif(($prot_name) && ($line =~ /^\w/)){

		my ($locus,$prev_annot,$seq_score,$seq_access,$fs_score,$fs_access,$fs_db,$fs_model,$gs_score,$gs_access,$gs_db,$gs_model) = split("\t",$line);
		
		if ($fs_score ne "-"){
			my $db_loc = $provided_struct{$fs_db};
			my $provided_pdb = $locus;
			unless ($fs_model eq "-"){
				$provided_pdb .= "-".$fs_model;
			}
			push(@{$results{$prot_name}},[$db_loc,$provided_pdb,$fs_access]);
		}
		elsif ($gs_score ne "-"){
			my $db_loc = $provided_struct{$gs_db};
			my $provided_pdb = $locus;
			unless ($gs_model eq "-"){
				$provided_pdb .= "-".$gs_model;
			}
			push(@{$results{$prot_name}},[$db_loc,$provided_pdb,$gs_access]);
		}

	}

}
close IN;

foreach my $protein (keys(%results)){
	print ("Working on visualizations for $protein...\n");
	foreach my $match ((@{$results{$protein}})){
		my ($db_loc,$provided_pdb,$match_pdb) = @{$match};
		my $temp_provided_pdb = "$outdir/$protein/$provided_pdb.temp.pdb";
		# print($temp_provided_file."\n");
		my $provided_pdbs = `ls $db_loc/$provided_pdb*.pdb*`;
		my @provided_pdbs = split("\n",$provided_pdbs);
		if ($provided_pdbs[0] =~ /\.gz$/){
			system ("zcat $db_loc/$provided_pdb* > $temp_provided_pdb");
		}
		else{
			system ("cat $db_loc/$provided_pdb* > $temp_provided_pdb");
		}

		my $temp_match_pdb = "$outdir/$protein/$match_pdb.temp.pdb";
		my $match_pdbs = `ls $uniprot_struct/$match_pdb*.pdb*`;
		my @match_pdbs = split("\n",$match_pdbs);
		if ($match_pdbs[0] =~ /\.gz$/){
			system ("zcat $uniprot_struct/$match_pdb*.pdb* > $temp_match_pdb");
		}
		else{
			system ("cat $uniprot_struct/$match_pdb*.pdb* > $temp_match_pdb");
		}

		my $cxs_name = "$outdir/$protein/${provided_pdb}_${match_pdb}.cxs";
		if (-e $cxs_name) { print "  Alignment between $provided_pdb and $match_pdb found. Skipping alignment...\n"; }
		else {
			# ChimeraX API calling
			print "  Aligning $provided_pdb to $match_pdb with ChimeraX\n";
			system (
				"chimerax 1>/dev/null --nogui $script \\
				  -p $temp_provided_pdb \\
				  -m $temp_match_pdb \\
				  -o $outdir/$protein\n
			") == 0 or checksig();
		}
	}

	system ("rm $outdir/$protein/*.temp.pdb");
}

### Subroutine(s)
sub checksig {

	my $exit_code = $?;
	my $modulo = $exit_code % 255;

	print "\nExit code = $exit_code; modulo = $modulo \n";

	if ($modulo == 2) {
		print "\nSIGINT detected: Ctrl+C => exiting...\n";
		exit(2);
	}
	elsif ($modulo == 131) {
		print "\nSIGTERM detected: Ctrl+\\ => exiting...\n";
		exit(131);
	}

}