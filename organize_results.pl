#!/usr/bin/perl
## Pombert Lab 2022

my $name = 'organize_results.pl';
my $version = '1.5.5';
my $updated = '2022-07-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}

USAGE		${name} \\
		  -m Queri3D/UNIPROT_SCRAP_RESULTS/metadata.log \\
		  -p Queri3D/GESAMT_Results_Parsed/All_Parsed_Results.matches

OPTIONS
-m (--metadata)	metadata.log from uniprot_scraper.py
-st (--struct)	File containing parsed 3D homology results
-sq (--seqnc)	File containing parsed Seq homology results
-a (--annot)	Optional: tab-delimited file containing annotations for predicted structures
-o (--outdir)	Output directory [Default: RESULTS]
EXIT

die "\n".$usage."\n\n" unless @ARGV;

my $struct_file;
my $seqnc_file;
my $metadata_file;
my $annot_file;
my $outdir = "RESULTS";

GetOptions(
	'st|struct=s' => \$struct_file,
	'sq|seqnc=s' => \$seqnc_file,
	'm|metadata=s' => \$metadata_file,
	'a|annot=s' => \$annot_file,
	'o|outfile=s' => \$outdir,
);

unless (-d $outdir){
	make_path($outdir,{mode=>0755});
}

###################################################################################################
## Acquire and Organize metadata
###################################################################################################

### METADATA
## KEY => UNIPROT ACCESSION
## [0] PROT_NAME
## [1] ORG_NAME
## [2] FASTA_FILE
## [3] FEATURES
## [4] PDB_FILES
my %metadata;

### STRUCTURE LINK
## KEY => PDBCODE_CHAIN
## [0] PROT_NAME
## [1] UNIPROT_ACCESSION
my %struct_link;

## Track proteins
my %proteins;

my %features;

my $accession;
open META, "<", $metadata_file or die ("Unable to open file $metadata_file: $!\n");
while (my $line = <META>){
	chomp ($line);
	START:
	unless($line eq ''){
		if ($line =~ /^>(\w+)/){
			$accession = $1;
		}
		elsif ($line =~ /^\tPROTEIN_NAME/){
			$line = <META>;
			chomp ($line);
			my ($prot_name) = $line =~ /\t\t(.*)/;
			$metadata{$accession}[0] = $prot_name;
			$proteins{$prot_name} = 1;
		}
		elsif ($line =~ /^\tORGANISM_NAME/){
			$line = <META>;
			chomp ($line);
			my ($org_name) = $line =~ /\t\t(.*)/;
			$metadata{$accession}[1] = $org_name;
		}
		elsif ($line =~ /^\tFASTA/){
			$line = <META>;
			chomp ($line);
			my ($fasta) = $line =~ /\t\t(.*)/;
			$metadata{$accession}[2] = $fasta;
		}
		elsif ($line =~ /^\tFEATURES/){
			while (1){
				$line = <META>;
				chomp ($line);
				if ($line =~ /^\t\t(.*)/){
					push(@{$metadata{$accession}[3]},$1);
				}
				else{
					if (join(";",@{$metadata{$accession}[3]}) ne "N\\A"){
						$features{$metadata{$accession}[0]} = join(";",@{$metadata{$accession}[3]});
					}
					goto START;
				}
			}
		}
		elsif ($line =~ /^\tSTRUCTURES/){
			while (1){
				$line = <META>;
				if ($line){
					chomp ($line);
					if ($line =~ /^\t\t(\w+)\t(.?)\t/){
						my $pdb = $1;
						my $chain = $2;
						if ($chain =~ /[a-zA-Z]/){
							push(@{$metadata{$accession}[4]},$pdb."_".$chain);
							@{$struct_link{$pdb."_".$chain}} = ($metadata{$accession}[0],$accession);
						}
						else{
							push(@{$metadata{$accession}[4]},$pdb);
							@{$struct_link{$pdb}} = ($metadata{$accession}[0],$accession);
						}
					}
					else{
						goto START;
					}
					
				}
				else{
					last;
				}
			}
		}
	}
}
close META;

###################################################################################################
## Acquire annotations
###################################################################################################
my %annotations;

if ($annot_file){
	open IN, "<", $annot_file or die "Unable to access file $annot_file: $!\n";
	while (my $line = <IN>){
		chomp($line);
		my ($locus,$annot,$notes) = split("\t",$line);
		$annotations{$locus} = $annot;
	}
}

###################################################################################################
## Acquire and Organize Sequence Results
###################################################################################################

### Tracks all results
## KEY1 => PROTEIN_NAME
## KEY2 => LOCUS_TAG
## [0] SeqHom data
## [1] 3DHom data
## [2] Combination of Eval and Qscore to sort
my %all_results;
my %seq_results;

## Keeps track of how many different proteins a locus hits against
my %loci_count;

open IN, "<", $seqnc_file or die "Unable to access $seqnc_file: $!\n";

while (my $line = <IN>){
	chomp($line);
	unless ($line eq "" || $line =~ /^###/){
		if ($line =~ /^## (\w+)/){
			$accession = $1;
		}
		else{

			## Data is going to be sorted by protein name to reduce entries; get it from metadata using accession
			my $prot_name = $metadata{$accession}[0];
			
			## Convert line to array
			my @data = split("\t",$line);
			
			## Remove locus from data
			my $locus = shift(@data);

			unshift(@data,$accession);

			$loci_count{"seq"}{$locus}++;

			## Store SeqHom data under the protein_name and locus_tag
			@{$all_results{$prot_name}{$locus}[0]} = @data;
			@{$seq_results{$prot_name}{$locus}} = @data;

			## Add the eval to score keeper
			$all_results{$prot_name}{$locus}[2] = -$data[9];
		}
	}
}

close IN;

###################################################################################################
## Acquire and Organize Structure Results
###################################################################################################

my %struct_results;

open IN, "<", $struct_file or die "Unable to access $struct_file: $!\n";

## Track PDB code
my $struct;

while (my $line = <IN>){
	chomp($line);
	unless($line eq "" || $line =~ /^###/){
		if ($line =~ /^## (\w+)/){
			$struct = $1;
		}
		else{

			my ($prot_name, $accession) = @{$struct_link{$struct}};
			
			## Convert line to array
			my @data = split("\t",$line);
			
			my $locus = shift(@data);
			unshift(@data,$struct);
			unshift(@data,$accession);

			## Use only the best hit per locus
			if ($all_results{$prot_name}{$locus}){
				unless ($all_results{$prot_name}{$locus}[1]){
					@{$all_results{$prot_name}{$locus}[1]} = @data;
					@{$struct_results{$prot_name}{$locus}} = @data;
					$all_results{$prot_name}{$locus}[2] += $data[3];
					$loci_count{"stc"}{$locus} ++;
				}
			}
			else{
				@{$all_results{$prot_name}{$locus}[1]} = @data;
				$all_results{$prot_name}{$locus}[2] = $data[3];
				$loci_count{"stc"}{$locus} ++;
			}
		}
	}
}
close IN;

###################################################################################################
## Print sequence results
###################################################################################################

open OUT, ">", $outdir."/sequence_results.tsv" or die "Unable to access file $outdir/sequence_results.tsv: $!\n";
print OUT "### LOCUS\tANNOTATION\tACCESSION\tPIDENT\tLENGTH\tMISMATCH\tGAPOPEN\tQSTART\tQEND\tSSTART\tSEND\tEVAL\tBITSCORE\n\n";
foreach my $prot (sort(keys(%seq_results))){
	print OUT "## $prot\n";
	foreach my $locus (sort{$seq_results{$prot}{$b}[9] <=> $seq_results{$prot}{$a}[9]}(keys(%{$seq_results{$prot}}))){
		print OUT $locus;
		if ($annotations{$locus}){
			print OUT "\t".$annotations{$locus}."\t";
		}
		else{
			print OUT "\t-";
		}
		print OUT "".(join("\t",@{$seq_results{$prot}{$locus}}))."\n";
	}
	print OUT "\n";
}
close OUT;

###################################################################################################
## Print structure results
###################################################################################################

open OUT, ">", $outdir."/structure_results.tsv" or die "Unable to access file $outdir/structure_results.tsv: $!\n";
print OUT "### LOCUS\tANNOTATION\tACCESSION\tPDB\tSOURCE\tQSCORE\tR.M.S.D.\tSEQID\tNalign\tnRES\n\n";
foreach my $prot (sort(keys(%struct_results))){
	print OUT "## $prot\n";
	foreach my $locus (sort{$struct_results{$prot}{$b}[3] <=> $struct_results{$prot}{$a}[3]}(keys(%{$struct_results{$prot}}))){
		if ($annotations{$locus}){
			print OUT $locus."\t".$annotations{$locus}."\t";
		}
		else{
			print OUT "\t-";
		}
		print OUT (join("\t",@{$struct_results{$prot}{$locus}}))."\n";
	}
	print OUT "\n";
}
close OUT;

###################################################################################################
## Print compiled results
###################################################################################################

open OUT, ">", $outdir."/compiled_results.tsv" or die "Unable to access file $outdir/compiled_results.tsv: $!\n";
my %loci_record;
# print OUT "### LOCUS\tANNOTATION\tSEQ_HOM_EVAL\tACCESSION\t(SEQ_HOM_HIT/SEQ_HOM_MATCHES)\t";
# print OUT "STC_HOM_QSCORE\tSTC_PDB\tPDB_SOURCE\t(STC_HOM_HIT/STC_HOM_MATCHES)\n\n";

print OUT "### LOCUS\tANNOTATION\tSEQ_HOM_EVAL\t";
print OUT "STC_HOM_QSCORE\tSTC_PDB\tPDB_SOURCE\n\n";

foreach my $prot (sort(keys(%proteins))){

	if ($all_results{$prot}){
		
		print OUT "## $prot\t(3D Features:";
		
		if ($features{$prot}){
			print OUT "[$features{$prot}])\n";
		}
		else{
			print OUT "[N/A])\n";
		}

		foreach my $locus (sort{$all_results{$prot}{$b}[2] <=> $all_results{$prot}{$a}[2]}(keys(%{$all_results{$prot}}))){

			my $score = $all_results{$prot}{$locus}[2];
			
			print OUT "$locus";
			
			if ($annotations{$locus}){
				print OUT "\t$annotations{$locus}";
			}
			else{
				print OUT "\t-";
			}

			if ($all_results{$prot}{$locus}[0]){
				$loci_record{"seq"}{$locus}++;
				my @seq_data = @{$all_results{$prot}{$locus}[0]};
				my ($seq_accession, $seq_pident, $seq_length, $seq_mismatch) = @seq_data[0..3];
				my ($seq_gapopen, $seq_qstart, $seq_qend, $seq_sstart) = @seq_data[4..8];
				my ($seq_send, $seq_eval, $seq_bitscore) = @seq_data[8..10];
				my $seq_org = $metadata{$seq_accession}[1];
				# print OUT "\t".$seq_eval."\t".$seq_accession."\t(".$loci_record{"seq"}{$locus}."/".$loci_count{"seq"}{$locus}.")";
				print OUT "\t".$seq_eval."\t".$seq_accession;
			}
			else{
				# print OUT "\t-"x3;
				print OUT "\t-"x2;
			}

			if ($all_results{$prot}{$locus}[1]){
				$loci_record{"stc"}{$locus}++;
				my @stc_data = @{$all_results{$prot}{$locus}[1]};
				my ($stc_accession, $stc_pdb, $stc_source, $qscore) = @stc_data[0..3];
				my ($stc_rmsd, $stc_seqid, $nalign, $nres) = @stc_data[4..7];
				# print OUT "\t".$qscore."\t".$stc_pdb."\t".$stc_source."\t(".$loci_record{"stc"}{$locus}."/".$loci_count{"stc"}{$locus}.")";
				print OUT "\t".$qscore."\t".$stc_pdb."\t".$stc_source;
			}
			else{
				# print OUT "\t-"x2;
				print OUT "\t-";
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
}
close OUT;