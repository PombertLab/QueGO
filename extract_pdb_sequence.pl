#!/usr/bin/perl
## Pombert Lab

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

my $name = 'extract_pdb_sequence.pl';
my $version = '0.1';
my $updated = '2022-04-11';
my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Extracts the amino acid sequence from pdb files.

USAGE		${name} \\
		  -p 1be3.pdb \\
		  -o EXTRACTED_FASTAS

OPTIONS
-p (--pdb)	PDB files to extract
-o (--out)	Directory to store extracted FASTA files [Default: EXTRACTED_FASTAS]
EXIT

die "\n$usage\n" unless @ARGV;

my @pdbs;
my $outdir = 'EXTRACTED_FASTAS';

GetOptions(
	'p|pdb=s{1,}' => \@pdbs,
	'o|out=s' => \$outdir,
);

unless (-d $outdir){
	make_path($outdir,{mode=>0755});
}

my %AAs;

initialize();

foreach my $file (@pdbs){

	open IN, "<", $file or die("Unable to read from $file: $!\n");

	my $fasta = "";

	while (my $line = <IN>){
		chomp($line);
		if($line =~ /^ATOM.{9}CA\s{2}(\w{3})/){
			$fasta .= $AAs{$1};
		}
	}

	close IN;

	my @fasta = unpack("(A60)*",$fasta);
	
	my ($filename) = $file =~ /(\w+)(?:\.\w+)+$/;
	
	open OUT, ">", "$outdir/$filename.fasta" or die("Unable to write to $file: $!\n");
	
	print OUT ">$filename\n";
	
	while (my $line = shift(@fasta)){
		print OUT $line."\n";
	}
}

sub initialize {
	%AAs = ('ALA' => 'A',
			'ASX' => 'B',
			'CYS' => 'C',
			'ASP' => 'D',
			'GLU' => 'E',
			'PHE' => 'F',
			'GLY' => 'G',
			'HIS' => 'H',
			'ILE' => 'I',
			'LYS' => 'K',
			'LEU' => 'L',
			'MET' => 'M',
			'ASN' => 'N',
			'PRO' => 'P',
			'GLN' => 'Q',
			'ARG' => 'R',
			'SER' => 'S',
			'THR' => 'T',
			'VAL' => 'V',
			'TRP' => 'W',
			'TYR' => 'Y',
			'GLX' => 'Z'
	);
}