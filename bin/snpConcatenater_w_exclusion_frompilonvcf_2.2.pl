#!/usr/bin/perl -w

#author Maha Farhat 

#Reads in all the vcf files in a directory, exclusion BED file, ID failed file and takes option INDEL and REGION. If REGION opion provided need to provide START= and STOP= coordinates.The region option results in an MSA of the full region with SNPs/indels introduced, and not just SNP concatenation
#for each line in the vcf file it excludes those lines that have a reference position between the start and end positions of an exclusion BED file 
# exports this to a single multiple alignment fasta file with the file name as the first field
#does some hard filtering of variants based on quality, excludes low


##example command: perl snpConcatenater_w_exclusion_frompilonvcf.pl -e <file_w_coord_start_stop_to_exclude in BED format> -f <failed_ID_list> -v [INDEL|SNP] -l [REGION|WHOLE] -h heterothresh[default is 0.1] -q qualitythresh[default is 20] -c <startcoord-zerobased>-<endcoord-1based> -s <strand:pos or neg> -o output_filename.fasta> -t coords_filename.txt\n REGION is optional; IF REGION is provided start stop and strand need to both be defined\n"
##note support for concurrent REGIOn & INDEl option is buggy and a work-in-progress use with caution


use warnings;
use strict;
use Getopt::Std;

my @tempFiles;

my $lsFileCommand = 'ls -1 *.vcf > files.txt';

system($lsFileCommand);

my @fileListRaw =&ReadInFile('files.txt');
my @fileList;

my $noargStatement="example command: perl snpConcatenater_w_exclusion_frompilonvcf.pl -e <file_w_coord_start_stop_to_exclude in BED format> -f <failed_ID_list> -v [INDEL|SNP] -l [REGION|WHOLE] -h heterothresh[default is 0.1] -q qualitythresh[default is 20] -c <startcoord-zerobased>-<endcoord-1based> -s <strand:pos or neg> -o output_filename.fasta> -t coords_filename.txt\n REGION is optional, but if provided start stop and strand need to be provided\n note support for concurrent REGIOn & INDEl option is buggy and a work-in-progress use with caution \n";
if ($#ARGV < 1) { 
    print STDERR "$noargStatement\n";
    die;
}

my %opts;
getopts('efv:l:cso:t:hq', \%opts);

my $option3=$opts{'v'};
my $option4=$opts{'l'};

my $regionstart; my $regionend; my $strand;

if ($option4 =~ m/REGION/i ) {
        my $i= $opts{'c'};
        if ($i =~ m/(\d+)-(\d+)/) {
                $regionstart=$1;
                $regionend=$2;
        }
	$strand= $opts{'s'};
        if ($regionstart !~ /\d+/) {
                print STDERR "$noargStatement\n";
                die;
        }
}
$regionstart++;

my @Start_excludedcoord;
my @End_excludedcoord;

if ($opts{'e'}) {
 my @ex =&ReadInFile($opts{'e'});
 foreach my $line (@ex) {
        chomp $line;
        my @pieces=split /\t/,$line;
        push (@Start_excludedcoord, $pieces[1]);
        push (@End_excludedcoord,$pieces[2]);
 }
}

foreach my $line (@fileListRaw) {
    chomp $line;
    if (length($line) > 0) {
        $line =~ s/\*//g;
        push (@fileList, $line);
    }
}


my %badIDs;
if ($opts{'f'}) {
 my @failed=&ReadInFile($opts{'f'});
 foreach my $line (@failed) {
    chomp $line;
    my @pieces=split /\s+/,$line;
    $badIDs{$pieces[0]}=1;
 }
}

push(@tempFiles, 'files.txt');

my %masterHash;
my %referenceHash;
my @allStrains;
my @allCoords;
my $qualThresh=$opts{'q'}||20;
my $heteroThresh=$opts{'h'}||0.10;
my @RvDNA;

if ($option4 =~ m/REGION/i) {
       my @sequence;
       @sequence = `/n/data1/hms/dbmi/farhat/bin/work-horse/bin/get_seq_coord.pl -coord ${regionstart}-${regionend} -nodefline /n/data1/hms/dbmi/farhat/bin/work-horse/bin/h37rv.fasta`;
       my $sequence;

       foreach my $line (@sequence) {
            chomp ($line);
            $sequence=$sequence.$line;
       }

       @RvDNA = split( '', $sequence);
       push (@Start_excludedcoord, 1);
       push (@End_excludedcoord, $regionstart-1 );
       my $i=$regionend+1;
       #print STDERR "$i\n";
       push (@Start_excludedcoord, $i);
       push (@End_excludedcoord, 4411532);
       
}


open(my $OUT1, '>', $opts{'t'});
open(my $OUT, '>', $opts{'o'});

FILE: for (my $i=0; $i<=$#fileList; $i++) {

    my @fileRaw =&ReadInFile("$fileList[$i]"); #("$fileList[$i]".'.out');

    my $strainName = $fileList[$i];
    $strainName =~ s/\.vcf//g;

    if (defined $badIDs{$strainName}) {
		print STDERR "skipping failed $strainName\n";
		next FILE;
    }
    push (@allStrains, $strainName);
    print STDERR "processing $strainName...\n";

    VAR: foreach my $line (@fileRaw) {
       chomp $line;
       if ($line =~ /^NC_000962\.3/) {

        	my @elements = split /\t/, $line;
  	        my ($from,$ref_allele,$allele)=($elements[1],@elements[3..4]);
		my $ex = &qualityControl($line);

		next VAR if $ex>0;

		if ( length($ref_allele) > 1 || length($allele) > 1) {
			if ($option3 =~ m/INDEL/i) {
	                 if (length($ref_allele) == length($allele)) { #multi-base substitution
					for (my $i=0; $i<length($ref_allele); $i++){
						my $pos=$from+$i;
						my $ref=substr($ref_allele,$i,1);
						my $base=substr($allele,$i,1);
						$masterHash{$strainName}{$pos} = $base;
						$referenceHash{$pos} = $ref;
				        push(@allCoords, $pos);
					}
                 } else {#indel
					if (length($ref_allele) > length($allele)) { #deletion
						for (my $i=0; $i<length($ref_allele); $i++){
							my $pos=$from+$i; my $base;
	                        my $ref=substr($ref_allele,$i,1);
							if (length($allele) <$i+1) {
								$base='-';
							} else {
        	                   	$base=substr($allele,$i,1);
							}
                	        $masterHash{$strainName}{$pos} = $base;
                        	$referenceHash{$pos} = $ref;
                            push(@allCoords, $pos);
						}
					} else { #insertion
						for (my $i=0; $i<length($allele); $i++){
                                                        my $pos=4411532+$from+$i;
							#if (defined $insSEENbefore{$pos}) {
							#	$pos=$pos+1000000;
							#}
							#$insSEENbefore{$pos}=1;
							my $base=substr($allele, $i, 1);
							if (length($ref_allele) <$i+1) {
								$referenceHash{$pos} ='-';
								$masterHash{$strainName}{$pos} = $base;
								push(@allCoords, $pos);
							} else {
								my $h=$from+$i;
								my $ref=substr($ref_allele,$i,1);
								$masterHash{$strainName}{$h} = $base;
								$referenceHash{$h} = $ref;
								push(@allCoords, $h);
							}
                        }
					}
		 		}
			} else {
                #print STDERR "excluded non single SUB\n";
	        	next VAR;
			}
        } else {
		#print STDERR "strain is $strainName ref is $ref_allele var is $allele \n";
    		$masterHash{$strainName}{$from} = $allele;
		#print STDERR "$masterHash{$strainName}{$from}\n";
       		$referenceHash{$from} = $ref_allele;
        	push(@allCoords, $from);
        }
       }
    }
}

my @uniqueCoords = &duplicateEliminator(@allCoords);
@uniqueCoords = sort {$a <=> $b} @uniqueCoords;

unless (exists($uniqueCoords[0])) {
	die "Warning: no variants in region found! Note: No fasta file written\n";
}

$"="\n";
print $OUT1 "Variants found at H37Rv coodinates:\n@uniqueCoords\n";
close $OUT1;
$"="";

if ($option4 =~ m/REGION/i) { #expand uniqueCoords to include all sites of the region, and place insertions in right place i.e. not at the end as above code does
	
	#print STDERR "original gene sequence is @RvDNA\n";
	foreach my $pos (@uniqueCoords) { #this loop adds gaps where any insertion has occurred relative to H37Rv such there reference RvDNA contains all the gaps
		if ($pos > 4411532) { 
			my $h= $pos - ( 4411532 + $regionstart + 1); #subtracting additional 1 as the dna sequence is stored in a zero based perl array 
			@RvDNA = (@RvDNA[0..$h],'-',@RvDNA[($h+1)..$#RvDNA]);
		}
	}
	#print STDERR "after adding gaps this becomes @RvDNA\n";
	
	for (my $i=0; $i<=$#allStrains; $i++) {
	    my $m=$allStrains[$i];
	    (my $y = $m) =~ s/-/_/g;
	    if ($y=~ /^\d+/) {
			print $OUT ">$y\n";
	    } else {
			print $OUT ">$y\n";
		}
		my @sequence=@RvDNA; 
		foreach my $pos (@uniqueCoords) {
			my $p;
			if ($pos < 4411532) { #subsitution or deletion
				$p=$pos - $regionstart;
			} else {
				$p= $pos - (4411532 + $regionstart); 
			}
			if (exists $masterHash{$m}{$pos} ) {
				#print STDERR "strain $m has a variant at relative position $p\t absolute position $pos\t reference base is $referenceHash{$pos} $sequence[$p] and the base in this strain is  $masterHash{$m}{$pos}\n";
				$sequence[$p]=$masterHash{$m}{$pos};
				#print STDERR "the sequence becomes: $sequence[$p]\n";
			}
		}
		if ($strand =~ m/neg/i) {
			@sequence = &revComp(@sequence);
		}
		$"="";
		print $OUT "@sequence\n";
	}
	print $OUT ">MT_H37Rv\n";
	if ($strand =~ m/neg/i) {
			@RvDNA = &revComp(@RvDNA);
	}
	print $OUT "@RvDNA\n";

} else { #ALL GENOME VARIANTS INCLUDED ONLY SNPCONCATENATE ALIGNMENTS ARE PRODUCED

	
	for (my $i=0; $i<=$#allStrains; $i++) {
	    my $m=$allStrains[$i];
	    $m=~ s/-/_/g;
	    if ($m=~ /^\d+/) {
			print $OUT ">$m\n";
	    } else {
			print $OUT ">$m\n";
	    }
	    for (my $j=0; $j<$#uniqueCoords; $j++) {
	        if (exists $masterHash{$allStrains[$i]}{$uniqueCoords[$j]}) {print $OUT $masterHash{$allStrains[$i]}{$uniqueCoords[$j]};}
	        else { print $OUT $referenceHash{$uniqueCoords[$j]}; }
	    }
	    print $OUT "\n";
	}

	print $OUT ">MT_H37Rv\n";
	for (my $j=0; $j<$#uniqueCoords; $j++) {print $OUT $referenceHash{$uniqueCoords[$j]};}
	print $OUT "\n";
}

close $OUT;
print STDERR "Done!\n";
foreach my $rmFile (@tempFiles) { system('rm '.$rmFile); }

########################################################################################
#SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES#
########################################################################################

#---------------------------------------------------------------------------------------
# Reads in a file and stores all the elements to an array
#---------------------------------------------------------------------------------------
sub ReadInFile
{
    my @FileName = @_;
    my @FileContents;
    open (FILE_OF_INTEREST, $FileName[0]) or die ("Unable to open the file called $FileName[0]");
    @FileContents = <FILE_OF_INTEREST>;
    close FILE_OF_INTEREST;
    return @FileContents;
}

#---------------------------------------------------------------------------------------
# Writes out a file
#---------------------------------------------------------------------------------------
sub WriteOutFile
{
  my @fileName = shift @_;
  my @fileContents =  @_;
  open (FILE_TO_WRITE_OUT, ">$fileName[0]") or die ("Unable to open the file called $fileName[0]");
  print FILE_TO_WRITE_OUT "@fileContents";
  close FILE_TO_WRITE_OUT;
}

#---------------------------------------------------------------------------------------
# Checks for duplications
#---------------------------------------------------------------------------------------
sub duplicateEliminator
{
    my @duplicateList = @_;
    my %uniqueHash;
    foreach my $i (@duplicateList) { $uniqueHash{$i} = 0;}
    my @uniqueList = keys(%uniqueHash);
    return @uniqueList;
}

#---------------------------------------------------------------------------------------
#  Check the qualtiy of a variant
#---------------------------------------------------------------------------------------

sub qualityControl
{
	my $line = shift @_;
	my $ex=0;
	my @elements = split /\t/, $line;
	my ($from,$ref_allele,$allele,$snpqual,$filter,$info)=($elements[1],@elements[3..7]);
	$ex=1 if $info =~ /IMPRECISE/i;
	(my $depth) = ($info =~ /DP=(\d+)/);
	$ex=2 if $ref_allele =~ /N/; #ambiguous reference
	$ex=2 if $allele =~ /N/; #ambiguous/imprecise change
	$ex=3 if $allele =~ m/,|</; #heterogenous allele or <dup>
    	(my $hqr)= ($info =~ /AF=((\d+(\.\d*)?)|(\.\d+))$/)||0.7666666;
	#if ($hqr eq ".") { #for some indels this may happen
	#	$hqr=0.7666666; #to patch this error
	#}
	if ($snpqual eq ".") {
		$snpqual=21;
	}
	if ($filter =~ m/;/ || $snpqual <=$qualThresh) { #|| $filter !~ m/PASS/ ) { #$hqr <=$heteroThresh ||
                $ex=4;
        }
	if ($filter !~ m/PASS/) {
                if ($filter =~ m/Amb/i && $hqr >$heteroThresh && $depth>9) {
                        $ex=0;
                } else {
                        $ex=5;
                }
        }
    	INTER: for (my $i=$#Start_excludedcoord; $i>=0; $i=$i-1) {
	#print STDERR "$Start_excludedcoord[$i] $End_excludedcoord[$i]\n";
    	if (($from >= $Start_excludedcoord[$i]) && ($from <= $End_excludedcoord[$i])) {
			$ex=6;
			last INTER;
        } 
    }
    return $ex;
}

#---------------------------------------------------------------------------------------
#  Check the qualtiy of a variant
#---------------------------------------------------------------------------------------
sub revComp 
{
	my $reversecomplement="";
	my @d= @_;
	foreach my $nucleotide(reverse(@d)) {
		if ($nucleotide =~ /a/i) {
			$reversecomplement.="T";
		} elsif ($nucleotide =~ /t/i) {
			$reversecomplement.="A";
		} elsif ($nucleotide =~ /g/i) {
			$reversecomplement.="C";
		} elsif ($nucleotide =~ /c/i) {
			$reversecomplement.="G";
		} elsif ($nucleotide =~ /-/) {
			$reversecomplement.="-";
		} else {
			die "$0:  Bad nucleotide!  [$nucleotide]\n";
		}
	}
	return split( '', $reversecomplement);
}
