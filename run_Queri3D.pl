#!/usr/bin/perl
## Pombert Lab 2022
my $name = "run_Queri3D.pl";
my $version = "0.5f";
my $updated = "2022-05-05";

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
SYNOPSIS	The purpose of this script is to download 3D structures related to given keywords
		from https://www.uniprot.org/ and runs homology searches of the known structures against
		an archive of predicted structures.

USAGE		$name \\
		  -g "telomere" \\
		  -v \\
		  -m "X-ray" \\
		  -p E_cuniculi_3D_structs \\
		  -t 4 \
		  -o Queri3D_telomere_Results

OPTIONS
-g (--go_annotation)		Search using gene ontolgy keyword
-v (--verified_only)		Search for genes that have been verified by UniProt
-m (--method)			Method used to obtain structure [Default = All] (i.e., X-ray, NMR)
-p (--predictions)		Directories containing predicted protein structures
-a (--archives)		GESAMT created archives
-t (--threads)			Number of threads to use [Default = 4]
-o (--outdir)			Output directory [Default = Queri3D_Results]
EXIT

die("\n$usage\n") unless(@ARGV);

my ($script,$pipeline_dir) = fileparse($0);
my $scraper_script = "$pipeline_dir/uniprot_scraper.py";
my $headers_script = "$pipeline_dir/PDB_headers.pl";
my $gesamt_script = "$pipeline_dir/run_GESAMT.pl";
my $parser_script = "$pipeline_dir/parse_GESAMT_results.pl";
my $metadata_script = "$pipeline_dir/add_metadata_to_results.pl";
my $custom;

my $go_annotation;
my $verified_only;
my $method;
my @predictions;
my @archives;
my $threads = 4;
my $outdir = "Queri3D_Results";

GetOptions(
	'g|go_annotations=s' => \$go_annotation,
	'v|verfied_only' => \$verified_only,
	'm|method=s' => \$method,
	'p|predictions=s{1,}' => \@predictions,
	'a|archive=s{1,}' => \@archives,
	't|threads=s' => \$threads,
	'o|outdir=s' => \$outdir,
	'c|custom=s' => \$custom, ## shhh, this is a secret tool for debugging purposes
);

my $gesamt_dir = "$outdir/GESAMT";
my $archives_dir = "$gesamt_dir/ARCHIVES";
my $results_dir = "$gesamt_dir/RESULTS";
my $uniprot_dir = "$outdir/UNIPROT_SCRAP_RESULTS";
my $pdb_dir = "$uniprot_dir/PDBs";
my $fasta_dir = "$uniprot_dir/FASTA";
my @dirs = ($outdir,$gesamt_dir,$archives_dir,$results_dir,$pdb_dir,$fasta_dir);

foreach my $dir (@dirs){
	unless(-d $dir){
		make_path($dir,{mode => 0755});
	}
}

my $start = localtime();
open LOG, ">", "$outdir/run_Queri3D.log" or die("Unable to open $outdir/run_Queri3D.log: $!\n");
print LOG "$0 started on $start\n";

$start = localtime();
print LOG "\nUniProt scraping started on $start\n";

unless(-f "$uniprot_dir/metadata.log"){

	my $flags = "";
	
	if($go_annotation){
		$flags .= "--go_annotation $go_annotation ";
	}


	if($method){
		$flags .= "--method $method ";
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

}
else{
	print "UniProt download completed previous! Moving to Homology Search!\n";
}

my $stop = localtime();

unless(-f "$uniprot_dir/metadata.log"){
	print LOG "\nUniProt scraping failed on $stop";
	exit;
}

print LOG "\nUniProt scraping completed on $stop\n";

if(@predictions||@archives){

	###################################################################################################
	## run_GESAMT.pl make archives
	###################################################################################################
	
	$start = localtime();
	print LOG "\nArchive creation started on $start\n";

	my %archives;

	foreach my $archive (@archives){
		my ($archive_path) = abs_path($archive);
		my ($archive_name) = $archive_path =~ /\/(\w+)\/*$/;
		unless(-d "$archives_dir/$archive_name"){
			system "ln -s $archive_path $archives_dir/$archive_name";
		}
		$archives{$archive_name} = "$archives_dir/$archive_name";
	}

	foreach my $predictor_dir (@predictions){
		my ($predictor) = $predictor_dir =~ /\/(\w+)\/*$/;
		unless(-d "$archives_dir/$predictor"){
			system "$gesamt_script \\
					-cpu $threads \\
					-make \\
					-arch $archives_dir/$predictor \\
					-pdb $predictor_dir
			";
		}
		else{
			print "Archive already created!\n";
		}
		$archives{$predictor} = "$archives_dir/$predictor";
	}


	$stop = localtime();
	print LOG "\nArchive creation completed on $stop\n";

	###################################################################################################
	## run_GESAMT.pl query search
	###################################################################################################

	$start = localtime();
	print LOG "\nGESAMT searches started on $start\n";

	my @results;

	foreach my $predictor (sort(keys(%archives))){
		my $archive = $archives{$predictor};
		unless(-d "$results_dir/$predictor"){
			opendir(QUERY,$pdb_dir);
			foreach my $file (readdir(QUERY)){
				unless(-d "$pdb_dir/$file"){
					my $pdb_name = fileparse($file,".pdb");
					unless(-f "$results_dir/$predictor/$pdb_name.normal.gesamt.gz"){
						system "$gesamt_script \\
								-cpu $threads \\
								-query \\
								-arch $archive \\
								-input $pdb_dir/$file \\
								-o $results_dir/$predictor \\
								-mode normal
						";
					}
				}
			}
		}
		push(@results,"$results_dir/$predictor");
	}

	$stop = localtime();
	print LOG "\nGESAMT searches completed on $stop\n";

	###################################################################################################
	## parse_GESAMT_results.pl query search
	###################################################################################################

	$stop = localtime();
	print LOG "\nStarted parsing results on $stop\n";

	system "$parser_script \\
			--results @results \\
			--outdir $results_dir
	";

	$stop = localtime();
	print LOG "\nFinished parsing results on $stop\n";

	###################################################################################################
	## add_metadata_to_results.pl query search
	###################################################################################################

	$start = localtime();
	print LOG "\nStarted adding metadata to results on $start\n";

	system "$metadata_script \\
			--metadata $uniprot_dir/metadata.log \\
			--parsed $results_dir/All_Parsed_Results.matches \\
			--outfile $results_dir/Parsed_Results_w_metadata.tsv
	";

	$stop = localtime();
	print LOG "\nFinished adding metadata to results on $stop\n";
}

$stop = localtime();
print LOG "\n$name completed on $stop\n";