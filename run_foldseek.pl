#!/usr/bin/perl
## Pombert Lab 2022
my $version = '0.2.3';
my $name = 'run_foldseek.pl';
my $updated = '2022-07-22';

use strict;
use warnings;
use File::Find;
use File::Basename;
use File::Path qw(make_path);
use POSIX 'strftime';
use Getopt::Long qw(GetOptions);

my @command = @ARGV; ## Keeping track of command line for log

## Usage definition
my $USAGE = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Run foldseek searches against PDB structures

REQUIREMENTS	Foldseek - https://github.com/steineggerlab/foldseek

CREATE DB	${name} -create -db /media/Data_2/FSEEK/rcsb -pdb /media/Data_2/PDB/
QUERY DB	${name} -query -db /media/Data_2/FSEEK/rcsb -input *.pdb -o ./FSEEK_RESULTS -z

OPTIONS:
-d (--db)	Foldseek database to create or query
-t (--threads)	CPU threads [Default: 12]
-v (--verbosity)	Verbosity: 0: quiet, 1: +errors, 2: +warnings, 3: +info [Default: 3]
-l (--log)	Log file [Default: foldseek.log]

## Creating a Foldseek database
-c (--create)	Create a foldseek database
-p (--pdb)	Folder containing the PDB files for the database

## Querying a Foldseek database
-q (--query)	Query a Foldseek database
-o (--outdir)	Output directory [Default: ./]
-i (--input)	PDF files to query
-a (--atype)	Alignment type [Default: 2]:
		  0: 3Di Gotoh-Smith-Waterman 
		  1: TMalign 
		  2: 3Di+AA Gotoh-Smith-Waterman
-m (--mseq)	Amount of prefilter sequences handed to the alignment [Default: 300]
-z (--gzip)	Compress output files [Default: off]

## Reference
van Kempen M, Kim S, Tumescheit C, Mirdita M, Söding J, and Steinegger M.
Foldseek: fast and accurate protein structure search. bioRxiv,
doi:10.1101/2022.02.07.479398 (2022)
OPTIONS
die "\n$USAGE\n" unless @ARGV;

## Defining options
my $db;
my $threads = 12;
my $verbosity = 3;
my $log = 'foldseek.log';

my $create;
my $pdb;

my $query;
my $outdir = './';
my @input;
my $atype = 2;
my $mseqs = 300;
my $gnuzip;
GetOptions(
	'd|db=s' => \$db,
	'l|log=s' => \$log,
	't|threads=i' => \$threads,
	'v|verbosity=i' => \$verbosity,
	'c|create' => \$create,
	'p|pdb=s' => \$pdb,
	'q|query' => \$query,
	'o|outdir=s' => \$outdir,
	'a|atype=i' => \$atype,
	'm|mseq=i' => \$mseqs,
	'i|input=s@{1,}' => \@input,
	'z|gzip' => \$gnuzip
);

## Creating log
my $date = strftime '%Y-%m-%d', localtime;
my $start = localtime();
my $tstart = time;
open LOG, ">>", "$log" or die "Can't create log file $log: $!\n";
print LOG "\nVERSION: $version\n"."COMMAND LINE: $name @command\n";
print LOG "Started on: $start\n";

## Program check
my $prog = `echo \$(command -v foldseek)`;
chomp $prog;
if ($prog eq ''){ 
	die "\nERROR: Cannot find foldseek. Please install foldseek in your path\n\n";
}

## Checking for unknown task
if (!defined $create and !defined $query){
	die "\nUnknown task. Please specify -create or -query on the command line.\n\n";
}

## Creating a foldseek database
if ($create){

	unless ($pdb){
		die "\nERROR: Please enter folder containing the PDB files for the database.";
	}

	my ($dbname,$dbpath) = fileparse($db);
	unless (-d $dbpath){ 
		make_path( $dbpath, { mode => 0755 } ) or die "Can't create folder $dbpath: $!\n";
	}

	system ("foldseek \\
			  createdb \\
			  --threads $threads \\
			  $pdb \\
			  $db") == 0 or checksig();
}

## Running foldseek queries/Skipping previously done searches
if ($query){

	unless (-d $outdir){ 
		make_path($outdir, {mode=>0755}) or die "Can't create folder $outdir: $!\n";
	}

	while (my $file = shift(@input)){

		my ($pdb, $dir) = fileparse($file);
		($pdb) = $pdb =~ /(\w+)\.pdb(?:\.gz)?$/;

		unless (-f "$outdir/$pdb.fseek" || -f "$outdir/$pdb.fseek.gz"){

			print "\n  Running foldseek on $file...\n";

			system ("foldseek \\
			  easy-search \\
			  --max-seqs $mseqs \\
			  --alignment-type $atype \\
			  --threads $threads \\
			  -v $verbosity \\
			  $file \\
			  $db \\
			  $outdir/$pdb.fseek \\
			  $outdir/tmp 1>/dev/null 2>$outdir/error.log") == 0 or checksig();
			
			if ($gnuzip){
				## Compressing data with GZIP to save some space
				if (-e "$outdir/$pdb.fseek"){
					system "gzip $outdir/$pdb.fseek";
				}
			}
		}

		## Searches can take a while, best to skip if done previously
		else { 
			print "  Skipping PDB file: $pdb => Foldseek result found in output directory $outdir\n";
		}
	}

	if ($gnuzip){
		print "\n  Results have been compressed with GZIP ...\n";
	}

	## Delete tmp directory
	if (-d "$outdir/tmp"){ system "rm -R $outdir/tmp"; }
}

my $end = localtime();
my $endtime = (time - $tstart)/60;
$endtime = sprintf ("%.2f", $endtime);
print LOG "Completed on: $end\n";
print LOG "Total run time: $endtime minutes\n";
close LOG;

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