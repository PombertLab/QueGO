#!/usr/bin/perl
## Pombert Lab 2022
my $name = "run_QueGO.pl";
my $version = "0.6.6";
my $updated = "2022-07-22";

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Basename;
use Cwd qw(abs_path);

my $usage = <<"EXIT";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	This script runs the QueGO pipeline from start to finish.

USAGE		$name \\
		  -k "telomere" \\
		  -v \\
		  -m "X-ray" "Predicted"\\
		  -s E_cuniculi_3D_structs \\
		  -t 4 \\
		  -o QueGO_telomere_Results

OPTIONS

## UNIPROT SCRAPER OPTIONS ##
-k (--go_keyword)	Search using gene ontolgy keyword
-v (--verified_only)	Search for genes that have been verified by experimental evidence
-m (--method)		Method used to obtain structure [Default = All] (X-ray, NMR, Predicted)
-u (--uniprot)		Previously performed UNIPROT_SCRAP_RESULTS

## HOMOLOGY OPTIONS ##
-s (--struct_sets)	Directories containing predicted protein structures
-h (--hom_tool)		3D Homology tool to use (FoldSeek or GESAMT) [Default: FoldSeek]
-a (--homology_arch)	3D homology archives (Archive must be compatible with --hom_tool)

## GENERAL OPTIONS ##
-t (--threads)		Number of threads to use [Default = 4]
-o (--outdir)		Output directory [Default = QueGO_Results]
EXIT

die("\n$usage\n") unless(@ARGV);

my $go_keyword;
my $verified_only;
my @method;
my @predictions;
my $hom_tool = "FOLDSEEK";
my @archives;
my $prot_fasta;
my $threads = 4;
my $uniprot;
my $outdir = "QueGO_Results";
my $custom;

GetOptions(
	'k|go_keyword=s' => \$go_keyword,
	'v|verfied_only' => \$verified_only,
	'm|method=s{1,}' => \@method,
	's|pred_struct=s{1,}' => \@predictions,
	'h|hom_tool=s' => \$hom_tool,
	'a|archive=s{1,}' => \@archives,
	'p|prot_fasta=s{1,}' => \$prot_fasta,
	't|threads=s' => \$threads,
	'u|uniprot=s' => \$uniprot,
	'o|outdir=s' => \$outdir,
	'c|custom=s' => \$custom, ## shhh, this is a secret tool for debugging purposes
);

## Setup script variables
my ($script,$pipeline_dir) = fileparse($0);
my $scraper_script = $pipeline_dir."/uniprot_scraper.py";
my $foldseek_script = $pipeline_dir."/run_foldseek.pl";
my $gesamt_script = $pipeline_dir."/run_GESAMT.pl";
my $parser_script = $pipeline_dir."/parse_3D_homology_results.pl";
my $metadata_script = $pipeline_dir."/organize_results.pl";

## Setup directory variables
my $uniprot_dir = $outdir."/UNIPROT_SCRAP_RESULTS";
my $fasta_dir = $uniprot_dir."/FASTA";
my $pdb_dir = $uniprot_dir."/PDBs";

my $seq_hom_dir = $outdir."/SEQUENCE_HOMOLOGY";

my $struct_hom_dir = $outdir."/STRUCTURE_HOMOLOGY";
my $arch_dir = $struct_hom_dir."/ARCHIVES";
my $struct_res_dir = $struct_hom_dir."/RESULTS";

my $results_dir = $outdir."/RESULTS";

my @dirs = (
	$outdir,

	$uniprot_dir,
	$fasta_dir,
	$pdb_dir,

	$seq_hom_dir,

	$struct_hom_dir,
	$arch_dir,
	$struct_res_dir,

	$results_dir
);

foreach my $dir (@dirs){
	unless(-d $dir){
		make_path($dir,{mode => 0755});
	}
}

###################################################################################################
## Getting UniProt data either from new WebScrap or old archive
###################################################################################################

unless (-f "$uniprot_dir/metadata.log"){
	### Use previously used scrap results
	if ($uniprot){
		if (-d "$uniprot"){
			system "cp -r $uniprot/* $uniprot_dir/";
		}
	}
	### Perform UniProt scraping
	elsif($go_keyword){
		### Check if scrap has been done previously
		
		my $flags = "";
		
		if($go_keyword){
			$flags .= "--go_keyword $go_keyword ";
		}


		if(@method){
			$flags .= "--method @method ";
		}

		if($verified_only){
			$flags .= "-v";
		}

		if($custom){
			$flags = "-c '$custom'";
		}

		system "$scraper_script \\
				--outdir $outdir/UNIPROT_SCRAP_RESULTS \\
				-df \\
				-ds \\
				-s \\
				$flags
		";

		unless(-f "$uniprot_dir/metadata.log"){
			die "\n[E]  UniProt scraping failed\n";
		}

	}
	else{
		die "\n[E]  Please provide a GO keyword or UNIPROT_SCRAP_RESULTS directory!\n";
	}
}
else{
	print "Previous UniProt scrap found at $uniprot_dir. Running on existing scrap results!\n";
}

my %archives;

###################################################################################################
## Copying precompiled 3D homology archives
###################################################################################################

if(@archives){
	
	### Copy pre-existing archives to new work enviroment
	foreach my $archive (@archives){

		opendir(ARCH,$archive);

		my $arch_tool = "FOLDSEEK";

		while (my $item = readdir(ARCH)){
			unless (-d $archive."/"."$item"){
				if ($item =~ /^gesamt.archive/){
					$arch_tool = "GESAMT";
					last;
				}
			}
		}

		close ARCH;

		if ($arch_tool eq $hom_tool){
			my ($archive_path) = abs_path($archive);
			my ($archive_name) = $archive_path =~ /\/(\w+)\/*$/;
			unless(-d $arch_dir."/".$hom_tool."/".$archive_name){
				system "cp -r $archive_path $arch_dir/$hom_tool/$archive_name";
			}
			$archives{$archive_name} = $archive_path;
		}
		else{
			continue;
		}

	}
}

###################################################################################################
## Creating 3D homology archives
###################################################################################################

ARCHIVE_CREATION:
if (@predictions){
	### Make FOLDSEEK archive from structure sets
	if (uc($hom_tool) eq "FOLDSEEK"){
		for my $structure_set (@predictions){
			my ($db_name) = $structure_set =~ /\/(\w+)\/*$/;
			unless (-d $arch_dir."/FOLDSEEK/".$db_name){	
				system ("
					$foldseek_script \\
					  --create \\
					  --db $arch_dir/FOLDSEEK/$db_name/$db_name \\
					  --pdb $structure_set \\
					  --threads $threads
				");
				$archives{$db_name} = $arch_dir."/FOLDSEEK/".$db_name;
			}
			else{
				continue;
			}
		}
	}
	### Make GESAMT archive from structure sets
	elsif (uc($hom_tool) eq "GESAMT"){
		foreach my $structure_set (@predictions){
			my ($db_name) = $structure_set =~ /\/(\w+)\/*$/;
			unless (-d $arch_dir."/GESAMT/".$db_name){
				system ("
					$gesamt_script \\
					  -cpu $threads \\
					  -make \\
					  -arch $arch_dir/GESAMT/$db_name \\
					  -pdb $structure_set
				");
				$archives{$db_name} =  $arch_dir."/GESAMT/".$db_name;
			}
			else{
				continue;
			}
		}
	}
	else{
		ASK:
		print "\n\n[W] $hom_tool is not a valid homology tool!\n\n Please select either GESAMT or FoldSeek:\t";
		chomp(my $selection = <STDIN>);
		if (uc($selection) eq "GESAMT" || uc($selection) eq "FOLDSEEK"){
			$hom_tool = $selection;
			goto ARCHIVE_CREATION;
		}
		else{
			goto ASK;
		}
	}
}

###################################################################################################
## Perform 3D homology searches
###################################################################################################

my @results;

for my $arch (keys(%archives)){
	my $arch_path = $archives{$arch};
	if (uc($hom_tool) eq "FOLDSEEK"){
		system ("
			$foldseek_script \\
			  --query \\
			  --db $arch_path/$arch \\
			  --input $pdb_dir/*\.pdb\.* \\
			  --outdir $struct_res_dir/FOLDSEEK/$arch \\
			  --gzip
		");
	}
	elsif (uc($hom_tool) eq "GESAMT"){
		system ("
			$gesamt_script \\
				--cpu $threads \\
				--query \\
				--arch $arch \\
				--input $pdb_dir/*\.pdb\.* \\
				--outdir $struct_res_dir/GESAMT/$arch \\
				--gzip \\
				-mode normal
		");
	}
}

###################################################################################################
## Parse 3D Homology Results
###################################################################################################

system "$parser_script \\
		--gesamt $struct_res_dir/GESAMT/* \\
		--foldseek $struct_res_dir/FOLDSEEK/* \\
		--outdir $results_dir
";

###################################################################################################
## Add metadata to GESAMT results
###################################################################################################

system "$metadata_script \\
		  --metadata $uniprot_dir/metadata.log \\
		  --foldseek $results_dir/FoldSeek_parsed_results.matches \\
		  --gesamt $results_dir/GESAMT_parsed_results.matches \\
		  --outfile $results_dir/compiled_results_w_metadata.tsv
";