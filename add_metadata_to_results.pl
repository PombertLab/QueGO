#!/usr/bin/perl
## Pombert Lab 2022

my $name = 'add_metadata_to_results.pl';
my $version = '0.5.1';
my $updated = '2022-05-24';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}

USAGE		${name} \\
		  -m Queri3D/UNIPROT_SCRAP_RESULTS/metadata.log \\
		  -p Queri3D/GESAMT_Results_Parsed/All_Parsed_Results.matches

OPTIONS
-m (--metadata)	metadata.log from uniprot_scraper.py
-p (--parsed)	File containing parsed 3D homology results
-a (--annot)	Optional: tab-delimited file containing annotations for predicted structures
-o (--outfile)	Output filename [Default: RESULTS_w_metadata.tsv]
EXIT

my $parsed_file;
my $metadata_file;
my $output = "results.tsv";

GetOptions(
	'p|parsed=s' => \$parsed_file,
	'm|metadata=s' => \$metadata_file,
	'o|outfile=s' => \$output,
);

my %metadata;
my %link;
my $accession;
my $organism;
my $prot_name;

open IN, "<", $metadata_file or die ("Unable to open file $metadata_file: $!\n");
while (my $line = <IN>){
	chomp ($line);
	unless($line eq ''){
		if ($line =~ /^\w+/){
			($accession,$organism,$prot_name) = split("\t",$line);
		}
		elsif ($line =~ /\t\t(\w+)\.pdb/){
			@{$metadata{$accession}} = ($accession,$organism,$prot_name);
			$link{$1} = $accession;
		}
	}
}
close IN;

my %results;
my $pdb_file;
open IN, "<", $parsed_file or die ("Unable to open file $parsed_file: $!\n");
while (my $line = <IN>){
	chomp ($line);
	if ($line =~ /\S/){
		if ($line =~ /^## (\w+)/){
			$pdb_file = $1;
			if ($link{$pdb_file}){
				$accession = $link{$pdb_file};
			}
			else{
				$accession = $pdb_file;
			}
		}
		elsif ($line =~ /^\w/) {
			my ($locus,@data) = split("\t",$line);
			if ($results{$accession}{$locus}){
				if ($results{$accession}{$locus}[2] > $data[1]){
					@{$results{$accession}{$locus}} = ($pdb_file,@data);
				}
			}
			else{
				@{$results{$accession}{$locus}} = ($pdb_file,@data);
			}
		}
	}
}

open OUT, ">", $output or die ("Unable to write to $output: $!\n");

foreach my $result (sort(keys(%results))){
	print OUT "### ".@{$metadata{$result}}[-1]."\t(Accession: ".@{$metadata{$result}}[0]."; ".@{$metadata{$result}}[1].")"."\n";
	foreach my $locus (sort{$results{$result}{$b}[2] <=> $results{$result}{$a}[2]}(keys(%{$results{$result}}))){
		print OUT $locus."\t".join("\t",@{$results{$result}{$locus}})."\n";
	}
	print OUT "\n"
}

close OUT;