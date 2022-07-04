#!/usr/bin/perl
## Pombert Lab, 2022

my $name = 'perform_sequence_search.pl';
my $version = '0.0.1';
my $updated = '2022-07-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

my $usage = <<"EXIT";


OPTIONS
-f (--faa)		.faa files from genome
-u (--uni)		UniProt FASTA directory
-t (--threads)	Threads [Default: 4]
-e (--eval)		E-value cutoff [Default: 1e-10]
-o (--outdir)	Output directory [Default: SEQUENCE_SEARCH]

EXIT

die "\n".$usage."\n\n" unless @ARGV;

my @subs;
my $queries;
my $threads = 4;
my $eval = "1e-10";
my $outdir = "SEQUENCE_SEARCH";

GetOptions(
	'f|fasta=s{1,}' => \@subs,
	'u|uni=s' => \$queries,
	't|threads=s' => \$threads,
	'e|eval=s' => \$eval,
	'o|out=s' => \$outdir,
);

my $results_dir = $outdir."/"."RESULTS";

my @dirs = ($outdir,$results_dir);

foreach my $dir (@dirs){
	unless(-d $dir){
		make_path($dir,{mode=>0755}) or die "Unable to create directory $dir: $!\n";
	}
}

unless (-f $outdir."/DB.dmnd"){
	system("diamond makedb --in @subs --db $outdir/DB");
}

opendir(DIR,$queries) or die "Unable to access directory $queries: $!\n";
foreach my $item (readdir(DIR)){
	unless (-d $queries."/".$item){
		my ($accession) = $item =~ /(\w+)\.fasta$/;
		print $accession."\n";
		unless (-f $results_dir."/".$accession.".diamond.6"){
			system("
				diamond \\
				blastp \\
				--threads $threads \\
				--db $outdir/DB \\
				--out $results_dir/$accession.diamond.6 \\
				--outfmt 6\\
				--query $queries/$item \\
				--evalue $eval \\
				1>/dev/null 2>$outdir/diamond.errors
			")
		}
	}
}

open OUT, ">", $outdir."/All_sequence_results.tsv" or die "Unable to access file $outdir/All_sequence_results.tsv: $!\n";
opendir(DIR,$results_dir) or die "Unable to access directory $results_dir: $!\n";
foreach my $item (readdir(DIR)){
	unless(-d $results_dir."/".$item){
		my ($accession) = $item =~ /(\w+).diamond.6/;
		open IN, "<", $results_dir."/".$item or die "Unable to access file $results_dir/$item: $!\n";
		print OUT "## $accession\n";
		while (my $line = <IN>){
			chomp($line);
			my @data = split("\t",$line);
			shift(@data);
			print OUT join("\t",@data)."\n";
		}
		print OUT "\n";
	}
}