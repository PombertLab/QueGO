#!/usr/bin/perl
## Pombert Lab 2022
my $name = "run_QueGO.pl";
my $version = "0.8.1";
my $updated = "2022-07-24";

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Basename;
use Cwd qw(abs_path);
use Term::ANSIColor;

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

## UNIPROT SCRAPER OPTIONS ##
-k (--go_keyword)	Search using gene ontolgy keyword
-v (--verified_only)	Search for genes that have been verified by experimental evidence
-n (--need_3D)	Require genes to have 3D structure
-m (--method)		Method used to obtain structure [Default = All] (X-ray, NMR, Predicted)
-u (--uniprot)		Previously performed UNIPROT_SCRAP_RESULTS

## SEQUENCE HOMOLOGY OPTIONS ##
-f (--fastas)		Files containing protein sequences (FASTAs extracted automatically from provided predicted structures if ignored)
-e (--eval)		E-value cut-off [Default: 1e-10]

## 3D HOMOLOGY OPTIONS ##
-s (--struct_sets)	Directories containing predicted protein structures (FASTAs extracted automatically if --fastas not provided)
-h (--hom_tool)		3D Homology tool to use (FoldSeek or GESAMT) [Default: FoldSeek]
-r (--homology_arch)	3D homology archives (Archive must be compatible with --hom_tool)
-t (--tmscore)		TM-score cut-off for FoldSeek [Default: 0.3]
-q (--qscore)		Q-score cut-off for GESAMT [Default: 0.3]

## GENERAL OPTIONS ##
-a (--annot)		TSV file containing existing annotations for predicted proteins
-w (--threads)		Number of threads to use [Default = 4]
-o (--outdir)		Output directory [Default = QueGO_Results]
EXIT

die("\n$usage\n") unless(@ARGV);

my @arguments;

for my $arg (@ARGV){
	push(@arguments,$arg);
}

my $master_start = time();
my $start;
my $stop;

my $go_keyword;
my $verified_only;
my $need_3D;
my @method;
my $uniprot;

my @prot_fasta;
my $seq_eval = 1e-10;

my @predictions;
my $hom_tool = "FOLDSEEK";
my @archives;
my $fs_tm = 0.3;
my $qscore = 0.3;

my $annot_file;
my $threads = 4;
my $outdir = "QueGO_Results";

my $custom;

GetOptions(
	'k|go_keyword=s' => \$go_keyword,
	'v|verfied_only' => \$verified_only,
	'n|need_3D' => \$need_3D,
	'm|method=s{1,}' => \@method,
	'u|uniprot=s' => \$uniprot,

	'f|fastas=s{1,}' => \@prot_fasta,
	'e|eval=s' => \$seq_eval,

	's|pred_struct=s{1,}' => \@predictions,
	'h|hom_tool=s' => \$hom_tool,
	'r|homology_arch=s{1,}' => \@archives,
	't|tmscore=s' => \$fs_tm,
	'q|qscore=s' => \$qscore,

	'a|annot=s' => \$annot_file,
	'w|threads=s' => \$threads,
	'o|outdir=s' => \$outdir,

	'c|custom=s' => \$custom, ## shhh, this is a secret tool for debugging purposes
);

## Setup script variables
my ($script,$pipeline_dir) = fileparse($0);
my $scraper_script = $pipeline_dir."/uniprot_scraper.py";
my $extract_script = $pipeline_dir."/extract_pdb_sequence.pl";
my $seq_hom_script = $pipeline_dir."/perform_sequence_search.pl";
my $foldseek_script = $pipeline_dir."/run_foldseek.pl";
my $mican_script = $pipeline_dir."/run_MICAN.pl";
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

print "\n\nSetting up working enviroment\n\n";

foreach my $dir (@dirs){
	unless(-d $dir){
		print "\tMaking directory $dir...\n";
		make_path($dir,{mode => 0755});
	}
	else{
		print "\tDirectory $dir already exists...\n";
	}
}

###################################################################################################
## Setting up log file
###################################################################################################

open LOG, ">", "$outdir/run_QueGO.log";
print LOG ($0);
for my $arg (@arguments){
	if(substr($arg,0,1) eq "-"){
		print LOG (" \\\n  $arg");
	}
	else{
		print LOG (" $arg")
	}
}
print LOG ("\n\n");

###################################################################################################
## Getting UniProt data either from new WebScrap or old archive
###################################################################################################

unless (-f "$uniprot_dir/metadata.log"){
	### Use previously used scrap results
	if ($uniprot){
		if (-d "$uniprot"){
			$start = time();
			print LOG "\n\tCopying UniProt scrap started at ".localtime($start)."\n";
			print "Utilizing previous UniProt scrap located at $uniprot...\n\n";
			system "cp -r $uniprot/* $uniprot_dir/";
			$stop = time();
			print LOG "\tUniProt scrap copying completed at ".localtime($stop)." (".duration($stop,$start).")\n";
		}
		else {
			print color 'red';
			print "\n\n[E]  Could not access previous UniProt scrap located at $uniprot...\n\n";
			print color 'reset';
			exit;
		}
	}
	### Perform UniProt scraping
	elsif($go_keyword){
		$start = time();
		print LOG "\n\tUniProt scrap started at ".localtime($start)."\n";
		print "\nStarting UniProt scrap...\n\n";
		my $flags = "";
		
		if($go_keyword){
			$flags .= "--go_keyword $go_keyword ";
		}

		if(@method){
			$flags .= "--method @method ";
		}

		if($need_3D){
			$flags .= "--structures ";
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
				$flags
		";

		unless(-f "$uniprot_dir/metadata.log"){
			print color 'red';
			print "\n\n[E]  UniProt scraping failed\n\n";
			print color 'reset';
			exit;
		}
		$stop = time();
		print LOG "\tUniProt scrap completed at ".localtime($stop)." (".duration($stop,$start).")\n";
	}
	else{
		print color 'red';
		print "\n\n[E]  Please provide a GO keyword or UNIPROT_SCRAP_RESULTS directory...\n\n";
		print color 'reset';
		exit;
	}
}
else{
	print color 'yellow';
	print "\nPrevious UniProt scrap found at $uniprot_dir. Utilizing previous results...\n";
	print color 'reset';
}

my %archives;

###################################################################################################
## Copying precompiled 3D homology archives
###################################################################################################

if(@archives){

	$start = time();

	print LOG "\n\tCopying archives started at ".localtime($start)."...\n";
	
	### Copy pre-existing archives to new work enviroment
	foreach my $archive (@archives){

		opendir(ARCH,$archive);

		my $arch_tool = "FOLDSEEK";

		## Checking if archive is GESAMT or FoldSeek
		while (my $item = readdir(ARCH)){
			unless (-d $archive."/"."$item"){
				if ($item =~ /^gesamt.archive/){
					$arch_tool = "GESAMT";
					last;
				}
			}
		}

		close ARCH;

		print("Provided archive at $archive is a $arch_tool archive...\n\n");

		if ($arch_tool eq $hom_tool){
			my ($archive_path) = abs_path($archive);
			my ($archive_name) = $archive_path =~ /\/(\w+)\/*$/;
			unless(-d $arch_dir."/".$hom_tool."/".$archive_name){
				system "cp -r $archive_path $arch_dir/$hom_tool/$archive_name";
			}
			else{
				print("\t$archive already exists in archive location $arch_dir/$hom_tool. Skipping...\n\n");
			}
			$archives{$archive_name} = $archive_path;
		}
		else{
			print color 'yellow';
			print("\t[W]  $archive is a $arch_tool archive and is not compatible with $hom_tool. Skipping...\n\n");
			print color 'reset';
		}

	}

	$stop = time();

	print LOG "\tArchive copying completed at ".localtime($start)." (".duration($stop,$start).")\n";
}

###################################################################################################
## Creating 3D homology archives
###################################################################################################

ARCHIVE_CREATION:
if (@predictions){

	$start = time();

	print LOG "\n\t$hom_tool archive creation started at ".localtime($start)."\n";
	print "\nCreating archives...\n";
	### Make FOLDSEEK archive from structure sets
	if (uc($hom_tool) eq "FOLDSEEK"){
		print "\tCreating archives for FoldSeek...\n";
		for my $structure_set (@predictions){
			my ($db_name) = $structure_set =~ /\/(\w+)\/*$/;
			$archives{$db_name} = $arch_dir."/FOLDSEEK/".$db_name;
			unless (-d $arch_dir."/FOLDSEEK/".$db_name){	
				print "\tCreating archive for $db_name...\n";
				system ("
					$foldseek_script \\
					  --create \\
					  --db $arch_dir/FOLDSEEK/$db_name/$db_name \\
					  --pdb $structure_set \\
					  --threads $threads
				");
			}
			else{
				print "\tArchive for $db_name already exists at $arch_dir/FOLDSEEK/...\n";
			}
		}
	}
	### Make GESAMT archive from structure sets
	elsif (uc($hom_tool) eq "GESAMT"){
		print "\tCreating archives for GESAMT...\n";
		foreach my $structure_set (@predictions){
			my ($db_name) = $structure_set =~ /\/(\w+)\/*$/;
			$archives{$db_name} =  $arch_dir."/GESAMT/".$db_name;
			unless (-d $arch_dir."/GESAMT/".$db_name){
				print "\tCreating archive for $db_name...\n";
				system ("
					$gesamt_script \\
					  -cpu $threads \\
					  -make \\
					  -arch $arch_dir/GESAMT/$db_name \\
					  -pdb $structure_set
				");
			}
			else{
				print "\tArchive for $db_name already exists at $arch_dir/GESAMT/...\n";
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

	$stop = time();

	print LOG "\t$hom_tool archive creation completed at ".localtime($stop)." (".duration($stop,$start).")\n";

}

###################################################################################################
## Extract PDB amino acid sequences
###################################################################################################
unless (@prot_fasta){
	$start = time();
	print LOG "\n\tProtein sequence extraction started at ".localtime($start)."\n";
	print "\nExtracting protein sequences from PDB files...\n";
	foreach my $structure_set (@predictions){
		system "
			$extract_script \\
			  --pdb $structure_set/*.pdb* \\
			  --out $seq_hom_dir/FASTA
		";
	}
	system "cat $seq_hom_dir/FASTA/*.faa > $seq_hom_dir/proteins.faa";
	$stop = time();
	print LOG "\tProtein sequence extraction completed at ".localtime($stop)." (".duration($stop,$start).")\n";
}
else{
	print LOG "\n\tProtein fastas provided. Skipping extraction...\n";
	foreach my $fasta (@prot_fasta){
		system "cat $fasta >> $seq_hom_dir/proteins.faa";
	}
}

###################################################################################################
## Perform sequence homology searches
###################################################################################################

unless (-f "$seq_hom_dir/All_sequence_results.tsv"){
	$start = time();
	print LOG "\n\tSequence homology searches started at ".localtime($start)."\n";
	print "\nPerforming sequence homology searches...\n";
	system "
		$seq_hom_script \\
		--faa $seq_hom_dir/proteins.faa \\
		--uni $fasta_dir \\
		--threads $threads \\
		--eval $seq_eval \\
		--outdir $seq_hom_dir
	";
	$stop = time();
	print LOG "\tSequence homology searches completed at ".localtime($stop)."( ".duration($stop,$start).")\n";
}
else{
	print LOG "\n\tSequence homology searches performed previously. Skipping search...\n";
}

###################################################################################################
## Perform 3D homology searches
###################################################################################################

my @results;

if (scalar(keys(%archives))>0){
	$start = time();
	print LOG "\n\t3D homology searches with $hom_tool started at ".localtime($start)."\n";
	for my $arch (keys(%archives)){
		print "\nPerforming 3D homology searches on $arch archive with ";
		my $arch_path = $archives{$arch};
		if (uc($hom_tool) eq "FOLDSEEK"){
			print "FoldSeek...\n";
			system ("
				$foldseek_script \\
				--query \\
				--db $arch_path/$arch \\
				--input $pdb_dir/*\.pdb* \\
				--outdir $struct_res_dir/FOLDSEEK/$arch \\
				--gzip
			");

		}
		elsif (uc($hom_tool) eq "GESAMT"){
			print "GESAMT...\n";
			system ("
				$gesamt_script \\
					--cpu $threads \\
					--query \\
					--arch $arch \\
					--input $pdb_dir/*\.pdb* \\
					--outdir $struct_res_dir/GESAMT/$arch \\
					--gzip \\
					-mode normal
			");
		}
	}
	$stop = time();
	print LOG "\t3D homology searches completed at ".localtime($stop)." (".duration($stop,$start).")\n";
}

if (uc($hom_tool) eq "FOLDSEEK"){
	$start = time();
	print LOG "\n\tTMscore calculation started at ".localtime($start)."\n";
	print "\nCalculating TMscores for FoldSeek results with MICAN...\n";
	system ("
		$mican_script \\
			--results_dir $struct_res_dir \\
			--uniprot_pdb $pdb_dir \\
			--predict_dir @predictions
	");
	$stop = time();
	print LOG "\tTMscore calculation completed at ".localtime($stop)." (".duration($stop,$start).")\n";
}

###################################################################################################
## Parse 3D Homology Results
###################################################################################################

$start = time();
print LOG "\n\tParsing 3D homology results started at ".localtime($start)."\n";
print "\nParsing 3D homology results...\n";
system "$parser_script \\
		--gesamt $struct_res_dir/GESAMT \\
		--foldseek $struct_res_dir/FOLDSEEK_w_MICAN \\
		--qscore $qscore \\
		--tm $fs_tm \\
		--outdir $results_dir
";
$stop = time();
print LOG "\tParsing 3D homology results completed at ".localtime($stop)."( ".duration($stop,$start).")\n";

###################################################################################################
## Compiling results and adding metadata
###################################################################################################

$start = time();
print LOG "\n\tCompiling results started at ".localtime($start)."...\n";
print "\nCompiling all evidences and adding metadata...\n";
system "$metadata_script \\
		  --metadata $uniprot_dir/metadata.log \\
		  --foldseek $results_dir/FoldSeek_parsed_results.matches \\
		  --gesamt $results_dir/GESAMT_parsed_results.matches \\
		  --seqnc $seq_hom_dir/All_sequence_results.tsv \\
		  --annot $annot_file \\
		  --outfile $results_dir
";
$stop = time();
print LOG ("\tResult compilation completed on ".localtime($stop)." (".duration($stop,$start).")\n");

###################################################################################################
## End of Script
###################################################################################################

print color 'reset';
my $master_stop = time();
print LOG ("\n$0 completed on ".localtime($master_stop)." (".duration($master_stop,$master_start).")\n");

###################################################################################################
## Subroutines

sub duration {
	my $elapsed = ($_[0] - $_[1]);
	my $days = int($elapsed/(24*60*60));
	my $hours = $elapsed%(24*60*60);
	$hours = int($hours/(60*60));
	my $mins = $elapsed%(60*60);
	$mins = int($mins/60);
	my $secs = $elapsed%60;
	$secs = int($secs);
	my $duration = "Duration: $days days, $hours hours, $mins mins, $secs secs";
	return $duration;
}