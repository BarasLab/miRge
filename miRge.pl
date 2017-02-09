# miRge -- an ultrafast, rational program for small RNA-Seq alignment and analysis
# Copyright (C) 2014-2015 Alex Baras <cmitch48@jhmi.edu>
# Copyright (C) 2015 Chris Mitchell <cmitch48@jhmi.edu>
# Copyright (C) 2014-2015 Marc Halushka <mhalush1@jhmi.edu>
# Copyright (C) 2014 Jason Myers, jason_myers@urmc.rochester.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;

use Pod::Usage;
use Getopt::Long;
use List::Util qw(max);
use Time::HiRes qw/ time sleep/;
use Cwd 'abs_path';
use HTML::Table;
use GD::Graph::bars;
use GD::Graph::hbars;
use File::Basename;
use File::Spec;
use Data::Dumper;
use String::Util qw(trim);


## Path to programs, by default miRge will use its own copies of the public tools bowtie and cutadapt within the miRge.seqUtils folder
my $trim5 = 0;
my $trim3 = 0;
my $miRgePath = abs_path($0);
my $version = "2.0";
local $ENV{PATH} = "$ENV{PATH}:".$miRgePath;
$miRgePath =~ s/\/[^\/]+\.pl/\//;
my $refPath = $miRgePath."miRge.seqLibs/";
my $trimBinary = File::Spec->catfile($miRgePath,"trim_file.py");

if(not -d $refPath){
	mkdir $refPath;
}

my $settings={};
my $help;
my @sampleFiles;
my @sampleFileNames;
my @sampleNames;
my $phred64 = '';
my $isomirDiff = '';
my $versionAsk = '';

GetOptions($settings,('help' => \$help,'version' => \$versionAsk,'outputHTML=s','trim5=s','trim3=s','outputDir=s','adapter=s','species=s','CPU=s','SampleNames=s','SampleFiles=s','isomirCutoff=s', 'bowtie=s', 'phred64' => \$phred64, 'diff-isomirs' => \$isomirDiff));
$trim5 = $$settings{trim5}||"0";
$trim3 = $$settings{trim3}||"0";
@sampleFiles = split(',', $$settings{'SampleFiles'});

for (my $i=0; $i<(@sampleFiles); $i++) {
    $sampleFiles[$i] = trim($sampleFiles[$i]); 
	unless(-e $sampleFiles[$i]) {
		die "$sampleFiles[$i] cannot be found, please check the paths of the sample files.";
	}
	$sampleFiles[$i] =~ m/\/?([^\/]+$)/;
	$sampleFileNames[$i] = $1;
}

if(exists($$settings{'SampleNames'})){
    @sampleNames = split(',', $$settings{'SampleNames'});
}
else{
    @sampleNames = @sampleFiles; 
}

my $adapter = $$settings{adapter}||"none";
if($adapter eq 'illumina'){
	$adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAG';
}
elsif($adapter eq 'ion'){
	$adapter = '+11';
}

my $speciesType = $$settings{species}||"none";
my $isomirCutoff = $$settings{isomirCutoff}||0.9;
my $numCPU = $$settings{CPU}||1;
my $bowtieBinary = $$settings{bowtie}||"bowtie";
my $bwtIndexDir = File::Spec->catdir($refPath,$speciesType);

my $mirnaBWT = File::Spec->catdir($bwtIndexDir,"mirna");
my $hairpinBWT = File::Spec->catdir($bwtIndexDir,"hairpin");
my $otherBWT = File::Spec->catdir($bwtIndexDir,"other");
my $mrnaBWT = File::Spec->catdir($bwtIndexDir,"mrna");

my $seqHash = {}; # key is nucleotide sequence, each record is a hash with
		  # 'annot' a boolean array where element [0] is overall matched status and each element [i] the matched status from each respective round of alignment
		  # 'quant' a integer array of the counts from each sample
my $mirHash = {}; # key is annot string from miRNA library records
my $readLengthHash = {}; # key is length
my $logHash = {}; # two primary keys annot and quant
my $graphHash = {};
my $tStamp = int(time);
my $outputPath = $$settings{outputDir}||'miRge.'.$tStamp;
my $outputHTML = $$settings{outputHTML}||"report.html";
system('mkdir -p '.$outputPath);
system('mkdir -p '.$outputPath.'/graphs');

my $t;
my $annotNames = ['exact miRNA', 'hairpin miRNA', 'non miRNA/mRNA RNA', 'mRNA', 'isomiR miRNA'];

pod2usage({ -message => "miRge version $version", -verbose => 0}) if( $versionAsk );
pod2usage( -verbose => 1) if( $help );
pod2usage( -verbose => 1) if( @sampleFiles<1 );
pod2usage({ -message => "You must provide an integral number of cores to use.", -verbose => 0}) if( $numCPU =~ m/\D/ );
pod2usage({ -message => "You must provide a species type.", -verbose => 0}) if( $speciesType eq "none" );


my $mirge_start_time = time;
print "\nChecking for bowtie and indices ...\n";
$t = time;
checkBowtie();
$t = getTimeDelta($t,time);
print "Bowtie and indices present ($t sec).\n";

print "\nStarting quantitation pipeline ...\n";
$t = time;
runQuantitationPipeline();
$t = getTimeDelta($t,time);
print "All samples completed ($t sec).\n";

print "\nStarting annotation pipeline ...\n";
$t = time;
runAnnotationPipeline();
$t = getTimeDelta($t,time);
print "All annotation cycles completed ($t sec).\n";

print "\nSummarizing and tabulating results ...\n";
$t = time;
summarize();
miRNAmerge();
filter();
generateGraphs();
writeHtmlReport();
writeDataToCSV();
$t = getTimeDelta($t,time);
print "Summary Complete ($t sec)\n";
$t = getTimeDelta($mirge_start_time,time);
print "miRge Completed ($t sec)\n";

############################# sub-routines #################################
sub calcEntropy{
	my $arrRef = shift;
	my @arr = @{$arrRef};
	my $sum = sumArray(\@arr);
	my $entropy = 0;
	for (my $i=0;$i<scalar(@arr);$i++){
		if($arr[$i] > 1){
			my $freq = $arr[$i]/$sum;
			$entropy = $entropy + -1*$freq*log($freq)/log(2);
		}
	}
	return $entropy;
}

sub sumArray{
	my $arrRef = shift;
	my $sum = 0;
	my @arr = @{$arrRef};
	for(my $i=0;$i<scalar(@arr);$i++){
		$sum = $sum + $arr[$i];
	}
	return $sum;
}

sub getTimeDelta{
	my $start = $_[0];
	my $end = $_[1];
	return sprintf('%.2f', $end-$start);
}

sub checkBowtieIndex {
	my $path = $_[0];
	my @indexes = glob("$path.*.ebwt");
	if(scalar(@indexes)){
		return $path;
	}
	my $index_type = $_[2];
	if($$settings{$index_type}){
		$path = $$settings{$index_type};
	}
	if(!(-e $path)){
		die "$path is unable to be found, please check the path.";
	}
	unless(-e $bwtIndexDir){
		mkdir($bwtIndexDir);
	}
	my $builder = $_[1];
	my($filename, $dirs, $suffix) = fileparse($path);
	if($suffix == 'fasta' || $suffix == 'fa' || $suffix == 'fna'){
		my $file = File::Spec->catdir($bwtIndexDir, $index_type);
		my @files = glob("$file.*.ebwt");
		if(scalar(@files)){
			return $file;
		}
		system("$builder $path $file");
		my $bowtieMap = File::Spec->catdir($bwtIndexDir, "index.map");
		my $bm;
		open $bm, ">>", $bowtieMap;
		print $bm "$path $file \n";
		if(-e $file){
			die "Error building $file.";
		}
		return $file;
		
	}
	return -e $path;
}

sub checkBowtie {
	my $bowtieIndex = "bowtie-build";
	# if(-e $bowtieBinary){
		# my($filename, $dirs, $suffix) = fileparse($bowtieBinary);
		# $bowtieIndex = File::Spec->catdir($dirs, "bowtie-build");
		# if(!(-e $bowtieIndex)){
			# die "Bowtie-build cannot be found at $bowtieIndex";
		# }
	# }
	# else{
		# die "Bowtie cannot be found at $bowtieBinary";
	# }
	$mirnaBWT = checkBowtieIndex($mirnaBWT, $bowtieIndex, 'mirna');
	$hairpinBWT = checkBowtieIndex($hairpinBWT, $bowtieIndex, 'hairpin');
	$otherBWT = checkBowtieIndex($otherBWT, $bowtieIndex, 'other');
	$mrnaBWT = checkBowtieIndex($mrnaBWT, $bowtieIndex, 'mrna');
}

sub runQuantitationPipeline {
	my $samplePrefix;
	my $cleanedReads;
	my $cleanedReads2;
	my $collapsedReads;
	my $i;
	
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
		print "Processing $sampleFiles[$i] ";
		$$logHash{'quantStats'}[$i]{'filename'} = $sampleFiles[$i];
		$$logHash{'quantStats'}[$i]{'name'} = $sampleNames[$i];
		$samplePrefix = basename($sampleFiles[$i]);
		$samplePrefix =~ s/\.fastq// ;
		$cleanedReads = $outputPath."/".$samplePrefix.".trim.fastq";
		$cleanedReads2 = $outputPath."/".$samplePrefix.".2.trim.fastq";
	
		trimRaw($sampleFiles[$i], $cleanedReads2, $i);
		print "cpuTime-trim:$$logHash{'quantStats'}[$i]{'cpuTime-trim'}, ";
        system("cutadapt -u $trim5 -u -$trim3 $cleanedReads2 > $cleanedReads") == 0 or die "Error: Couldn't trim first $trim5 and last $trim3 bases!\n";
		quantReads($cleanedReads, $i);
		print "cpuTime-uniq:$$logHash{'quantStats'}[$i]{'cpuTime-uniq'}\n";
		system("rm $cleanedReads $cleanedReads2");
	}
}

sub trimRaw {
	my $infile = $_[0];
	my $outfile = $_[1];
	my $logfile = "$outputPath/" . basename($infile) . ".log";
	my $sampleIndex = $_[2];
	my $fileText;
	my $fh;
		
	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-trim'} = time;
	my $command = "--adapter=$adapter --threads=$numCPU --infile=$infile --outfile=$outfile --logfile=$logfile";
	my $phred64 = `python $trimBinary $command` == 64;
	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-trim'} = getTimeDelta($$logHash{'quantStats'}[$sampleIndex]{'cpuTime-trim'}, time);
	
	open $fh, "<", "$logfile";
	$fileText = join('',<$fh>);
	close $fh;
	system("rm $logfile");
	$fileText =~ m/Starting reads: (\d+)/;
	$$logHash{'quantStats'}[$sampleIndex]{'totalReads'} = $1;
	$fileText =~ m/Processed reads: (\d+)/;
	$$logHash{'quantStats'}[$sampleIndex]{'trimmedReads'} = $1;
}

sub quantReads {
	my $infile = $_[0];
	my $sampleIndex = $_[1];
	my $fh;
	my $line;
	my $seqKey;
	
	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-uniq'} = time;
	open $fh, "<", $infile;
	while ($$line[0] = <$fh>) {
		$$line[1] = <$fh>;
		$$line[2] = <$fh>;
		$$line[3] = <$fh>;		
		chomp($$line[1]);
		
		$$seqHash{$$line[1]}{'quant'}[$sampleIndex] += 1;
	}
	close $fh;
	
	foreach $seqKey (keys %{$seqHash}) {
		if (!defined($$seqHash{$seqKey}{'annot'})) {
			$$seqHash{$seqKey}{'annot'}[0] = 0;
			$$seqHash{$seqKey}{'length'} = length($seqKey);
		}
		$$readLengthHash{$$seqHash{$seqKey}{'length'}}[$sampleIndex] += $$seqHash{$seqKey}{'quant'}[$sampleIndex];
	}

	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-uniq'} = getTimeDelta($$logHash{'quantStats'}[$sampleIndex]{'cpuTime-uniq'}, time);
}

sub runAnnotationPipeline {
	my $i;
	my $alignmentStatus;
	my $alignmentResult;	   

	my $lengthFilters = [-26,25,0,0,0];
	my $bwtCmd = "$bowtieBinary -a --best --strata --threads $numCPU";
# change to use --best --strata -a based on http://rnajournal.cshlp.org/content/22/8/1120.full
	if($phred64){
		$bwtCmd = "$bwtCmd --phred64-quals";
	}
	
	my $bwtCmdLines = [
		"$bwtCmd $mirnaBWT -n 0 -f $outputPath/SeqToAnnot.fasta 1>$outputPath/SeqToAnnot.sam 2>$outputPath/SeqToAnnot.log",
		"$bwtCmd $hairpinBWT -n 1 -f $outputPath/SeqToAnnot.fasta 1>$outputPath/SeqToAnnot.sam 2>$outputPath/SeqToAnnot.log",
		"$bwtCmd $otherBWT -n 1 -f $outputPath/SeqToAnnot.fasta 1>$outputPath/SeqToAnnot.sam 2>$outputPath/SeqToAnnot.log",
		"$bwtCmd $mrnaBWT -n 0 -f $outputPath/SeqToAnnot.fasta 1>$outputPath/SeqToAnnot.sam 2>$outputPath/SeqToAnnot.log",
		"$bwtCmd $mirnaBWT -l 15 -5 1 -3 2 -n 2 -f $outputPath/SeqToAnnot.fasta 1>$outputPath/SeqToAnnot.sam 2>$outputPath/SeqToAnnot.log"
	]; 
	# -- ALIGNMENT 1 -- length < 26, up to 0 mismatch to miRNA
	# -- ALIGNMENT 2 -- length > 25, up to 1 mismatch to hairpin
	# -- ALIGNMENT 3 -- any length, up to 1 mismatch other RNA
	# -- ALIGNMENT 4 -- any length, up to 0 mismatch to mRNA
	# -- ALIGNMENT 5 -- any length, up to 2 mismatches with special 5 vs 3 prime considerations to miRNA
	
	for ($i=0; $i<5; $i++) {
		# write the records with no matches/annotations and within length criteria (0=none,<0=smaller,>0=greater) to a fasta file to be used by the alignment/matching cmd tools	
		writeSeqToAnnot($$lengthFilters[$i]);
		
		# run alignment
		print "Starting Annotation-$$annotNames[$i] ";
		$$logHash{'annotStats'}[$i+1]{'cpuTime'} = time;
		$alignmentStatus = system($$bwtCmdLines[$i]);
		$$logHash{'annotStats'}[$i+1]{'cpuTime'} = time - $$logHash{'annotStats'}[$i+1]{'cpuTime'}; # in seconds
		$alignmentStatus = $alignmentStatus >> 8;
		print sprintf("cpuTime: %.2f\n", $$logHash{'annotStats'}[$i+1]{'cpuTime'});

		# parse alignment completed without error
		if ($alignmentStatus==0) {
			parseBowtieLog("$outputPath/SeqToAnnot.log",$i+1);
			$alignmentResult = parseAlignment("$outputPath/SeqToAnnot.sam");
			updateAnnotHash($alignmentResult,$i+1);		
		}
		else {
			print "Alignment exited with none-zero status.\n";
			exit;
		}
	
		# remove tmp files
		system("rm $outputPath/SeqToAnnot.fasta $outputPath/SeqToAnnot.sam $outputPath/SeqToAnnot.log");
	}
}

sub writeSeqToAnnot {
	my $lengthFilter = $_[0];
	my $filename = "$outputPath/SeqToAnnot.fasta";
	my $seqKey;
	my $fh;

	open $fh, ">", $filename;
	foreach $seqKey (keys %{$seqHash}) {
		if ($$seqHash{$seqKey}{'annot'}[0] == 0) {
			if ($lengthFilter < 0) {
				if ($$seqHash{$seqKey}{'length'}<-$lengthFilter) {
					print $fh ">$seqKey\n$seqKey\n";
				}
			} elsif ($lengthFilter > 0) {
				if ($$seqHash{$seqKey}{'length'}>$lengthFilter) {
					print $fh ">$seqKey\n$seqKey\n";
				}
			} else {
				print $fh ">$seqKey\n$seqKey\n";
			}
		}
	}
	close $fh;
}	

sub parseBowtieLog {
	my $file = $_[0];
	my $alignmentIndex = $_[1];
	my $fh;
	my $fileText;
	
	open $fh, "<", $file;
	$fileText = join('',<$fh>);
	close $fh;
	
	$fileText =~ m/reads processed: (\d+)/;
	$$logHash{'annotStats'}[$alignmentIndex]{'readsProcessed'} = $1;
	
	$fileText =~ m/reads with at least one reported alignment: (\d+)/;
	$$logHash{'annotStats'}[$alignmentIndex]{'readsAligned'} = $1;
}

sub parseAlignment {
	my $file = $_[0];
	my $alignmentResult = {};
	my $fh;
	my $line;

	open $fh, "<", $file;

	$line = <$fh>;
	while ($line =~ m/^@/) { $line = <$fh>;}

	while ($line) {
		$line = [split("\t",$line)];
		$$alignmentResult{$$line[0]}[0] = $$line[2];
		$$alignmentResult{$$line[0]}[1] = $$line[1];
		$line = <$fh>;
	}
	close $fh;

	return($alignmentResult);
}

sub updateAnnotHash {
	my $alignmentResult = $_[0];
	my $alignmentIndex = $_[1];
	my $seqKey;

	foreach $seqKey (keys %{$alignmentResult}) {
		if ($$alignmentResult{$seqKey}[0] ne '*') {
			$$seqHash{$seqKey}{'annot'}[0] = 1;
			$$seqHash{$seqKey}{'annot'}[$alignmentIndex] = $$alignmentResult{$seqKey}[0];
		} else {
			$$seqHash{$seqKey}{'annot'}[$alignmentIndex] = $$alignmentResult{$seqKey}[0];
		}
	}	
}

sub summarize {
	my $seqKey;
	my $mirKey;
	my $i;
	my $mirList = "$bowtieBinary-inspect -n $mirnaBWT";
		
	$mirList = [split("\n",`$mirList`)];	
	for ($i=0; $i<scalar(@{$mirList}); $i++) {
		$$mirHash{$$mirList[$i]}{'quant'} = [];
		$$mirHash{$$mirList[$i]}{'iscan'} = [];
	}
	
	foreach $seqKey (keys %{$seqHash}) {		
		for ($i=0; $i<scalar(@sampleFiles); $i++) {
			if (defined($$seqHash{$seqKey}{'quant'}[$i])) {
				$$logHash{'quantStats'}[$i]{'trimmedUniq'}++;
				if (defined($$seqHash{$seqKey}{'annot'}[1])|defined($$seqHash{$seqKey}{'annot'}[5])) {
					$$logHash{'quantStats'}[$i]{'mirnaReads'} += $$seqHash{$seqKey}{'quant'}[$i];
					#$$logHash{'quantStats'}[$i]{'mirnaUniq'}++;
					if (defined($$seqHash{$seqKey}{'annot'}[1])) {
						$$mirHash{$$seqHash{$seqKey}{'annot'}[1]}{'quant'}[$i] += $$seqHash{$seqKey}{'quant'}[$i];
						$$mirHash{$$seqHash{$seqKey}{'annot'}[1]}{'iscan'}[$i] += $$seqHash{$seqKey}{'quant'}[$i];
					} else {
						$$mirHash{$$seqHash{$seqKey}{'annot'}[5]}{'quant'}[$i] += $$seqHash{$seqKey}{'quant'}[$i];
					}
				}
				elsif (defined($$seqHash{$seqKey}{'annot'}[2])) {
					$$logHash{'quantStats'}[$i]{'hairpinReads'} += $$seqHash{$seqKey}{'quant'}[$i];
				}
				elsif (defined($$seqHash{$seqKey}{'annot'}[3])) {
					$$logHash{'quantStats'}[$i]{'ornaReads'} += $$seqHash{$seqKey}{'quant'}[$i];
				}
				elsif (defined($$seqHash{$seqKey}{'annot'}[4])) {
					$$logHash{'quantStats'}[$i]{'mrnaReads'} += $$seqHash{$seqKey}{'quant'}[$i];
				}
				else {
					$$logHash{'quantStats'}[$i]{'remReads'} += $$seqHash{$seqKey}{'quant'}[$i];
				}
			}
		}
	}
}

sub miRNAmerge {
	my $fh;
	my $line;
	my $i;
	my $j;
	
	my $mergeFile = File::Spec->catfile(File::Spec->catdir($refPath, $speciesType), 'merges.csv');
	
	if(-e $mergeFile){
		open $fh, "<", $mergeFile;
		while ($line = <$fh>) {
			chomp($line);
			$line = [split(',',$line)];
			for ($i=1; $i<scalar(@{$line}); $i++) {
				for ($j=0; $j<scalar(@sampleFiles); $j++) {
					if ($$mirHash{$$line[$i]}{'quant'}[$j]>0) {
						$$mirHash{$$line[0]}{'quant'}[$j] += $$mirHash{$$line[$i]}{'quant'}[$j];
						$$mirHash{$$line[0]}{'iscan'}[$j] += $$mirHash{$$line[$i]}{'iscan'}[$j];
					}
				}
				delete($$mirHash{$$line[$i]});
			}
		}
		close $fh;
	}
	else{
		print "Cannot find merges file, skipping merge step.\n";
	}
}

sub filter {
	my $mirKey;
	my $i;
	my $minCount;
	my $iscanFilter = 2;
	
	foreach $mirKey (keys %{$mirHash}) {
		for ($i=0; $i<scalar(@sampleFiles); $i++) {
			if ($$mirHash{$mirKey}{'iscan'}[$i]<$iscanFilter) {
				$$mirHash{$mirKey}{'quant'}[$i] = 0;
			}
		}		
	}
	
	foreach $mirKey (keys %{$mirHash}) {
		for ($i=0; $i<scalar(@sampleFiles); $i++) {
			if ($$mirHash{$mirKey}{'quant'}[$i]>0) {
				$$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'} += $$mirHash{$mirKey}{'quant'}[$i];
				$$logHash{'quantStats'}[$i]{'mirnaUniqFiltered'}++;
			}
		}
	}	
	
	for ($i=0; $i<scalar(@sampleFiles); $i++) {
		if($$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'} == 0){
			print STDERR "Warning: No miRNA reads were found in sample $sampleFiles[$i]. Please check your files and provided arguments.\n";
		}
	}
}

sub getValidFilename {
	my $fileIndex = shift;
	my $filename = shift;
	my $proposedPath = shift;
	my $fileKey = $fileIndex.$filename;
	if(defined($$graphHash{$fileKey})){
		return $$graphHash{$fileKey};
	};
	if(-e File::Spec->catdir($proposedPath, $filename)){
		my $fileindex = 1;
		$filename = $filename."_".$fileindex;
		while(-e File::Spec->catdir($proposedPath, $filename)){
			$fileindex = $fileindex + 1;
			$filename = substr($filename, 0, length($filename)-length($fileindex)).$fileindex;
		}
	};
	my $finalPath = File::Spec->catdir($proposedPath, $filename);
	$$graphHash{$fileKey} = $finalPath;
	return $finalPath;
}

sub generateGraphs {
	my $lengthKey;
	my $graph;
	my $graphData;
	my $countMax;
	my $lengthMax;
	my $i;
	my $j;
	my $fh;
	
	$lengthMax = max (keys %{$readLengthHash});
	for ($i=0; $i<scalar(@sampleFiles); $i++) {
		undef $graphData;
		for ($j=1; $j<=$lengthMax; $j++) {
			$$graphData[$j] = $$readLengthHash{$j}[$i];
		}
		
		$graph = GD::Graph::bars->new(600,300);
		$countMax = max @{$graphData};
		$graph->set(title => "$sampleNames[$i] (based on $$logHash{'quantStats'}[$i]{'trimmedReads'} reads)",
			    x_label => 'Read Length',
			    x_label_position => 0.5,
			    x_label_skip => 5,
			    y_label => 'Counts',
			    y_label_position => 0.5,
			    y_max_value => 10*(1+int($countMax/10)),
			    y_tick_number => 10,
			    borderclrs => undef,
			    dclrs => ['blue']);
		my $readPathName = $sampleFileNames[$i].'.readDistribution.png';
		my $readPath = getValidFilename($i, $readPathName, File::Spec->catdir($outputPath,'graphs'));
		open $fh, ">", $readPath;
		print $fh $graph->plot([[0..$lengthMax],$graphData])->png;
		close $fh;

		$graph = GD::Graph::hbars->new(600,300);
		$graph->set(title => "$sampleNames[$i] (based on $$logHash{'quantStats'}[$i]{'trimmedReads'} reads)",
			    y_label => 'Percentage',
			    y_label_position => 0.5,
				y_number_format => "%.2f",
			    borderclrs => undef,
			    dclrs => ['blue'],
			    show_values => 1,
				values_format => "%.3f");
		$graphData = [['miRNA','mRNA','other ncRNA','miRNA hairpin','unaligned'],
		[$$logHash{'quantStats'}[$i]{'mirnaReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'mrnaReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'ornaReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'hairpinReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'remReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'}]];
		my $alignPathName = $sampleFileNames[$i].'.readAlignments.png';
		my $alignPath = getValidFilename($i, $alignPathName, File::Spec->catdir($outputPath,'graphs'));
		open $fh, ">", $alignPath;
		print $fh $graph->plot($graphData)->png;
		close $fh;
	}
}

sub writeHtmlReport {
	my $filename = ($$settings{outputHTML}) ? $outputHTML : $outputPath . "/" . $outputHTML;
	my $fh;
	my $i;
	my $annotTable = new HTML::Table(0,0);
	my $quantTable = new HTML::Table(0,0);	
	my $fileTable = new HTML::Table(0,0);	
	
	$quantTable->setClass('tableBlue');
	$quantTable->setWidth(1000);
	$quantTable->addRow(('File name&#40;s&#41;','Total Input Reads','Trimmed Reads<br>(all / unique)','All miRNA Reads / Filtered miRNA Reads','Unique miRNAs','Hairpin Reads','Other Noncoding RNA Reads','mRNA Reads','Remaining Reads','Read Length and Composition Figures'));
	for ($i=0; $i<scalar(@sampleFiles); $i++) {
		my $readPathName = $sampleFileNames[$i].'.readDistribution.png';
		my $readPath = getValidFilename($i, $readPathName, File::Spec->catdir($outputPath,'graphs'));
		my $alignPathName = $sampleFileNames[$i].'.readAlignments.png';
		my $alignPath = getValidFilename($i, $alignPathName, File::Spec->catdir($outputPath,'graphs'));
		$quantTable->addRow(($$logHash{'quantStats'}[$i]{'name'},
				     $$logHash{'quantStats'}[$i]{'totalReads'},				 
				     '<table><tr></tr><tr><td>'.$$logHash{'quantStats'}[$i]{'trimmedReads'}.'&nbsp;/&nbsp;'.$$logHash{'quantStats'}[$i]{'trimmedUniq'}.'</td>'.'</tr></table>',
				     $$logHash{'quantStats'}[$i]{'mirnaReads'}.'&nbsp;/&nbsp;'.$$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'},
				     $$logHash{'quantStats'}[$i]{'mirnaUniqFiltered'},
				     $$logHash{'quantStats'}[$i]{'hairpinReads'},
				     $$logHash{'quantStats'}[$i]{'ornaReads'},
				     $$logHash{'quantStats'}[$i]{'mrnaReads'},
				     $$logHash{'quantStats'}[$i]{'remReads'},
				     #'<table><tr></tr><tr><td class="thumbnail1"><img src="../'.$readPath.'" width="100px" height="50px"><span><img src="../'.$readPath.'"></span></td><td class="thumbnail2"><img src="../'.$alignPath.'" width="100px" height="50px"><span><img src="../'.$alignPath.'"></span></td></tr></table>'));
				     '<table><tr></tr><tr><td class="thumbnail1"><img src="./graphs/'.$readPathName.'" width="100px" height="50px"><span><img src="./graphs/'.$readPathName.'"></span></td><td class="thumbnail2"><img src="./graphs/'.$alignPathName.'" width="100px" height="50px"><span><img src="./graphs/'.$alignPathName.'"></span></td></tr></table>'));
		}										 
					 
	$annotTable->setClass('tableBlue');
	$annotTable->setWidth(600);
	$annotTable->addRow(('Annotation-Round','# Unique Seqs','cpuTime(sec)'));
	$annotTable->addRow(('all sequences',scalar(keys %{$seqHash}),));	
	for ($i=1; $i<=5; $i++) {
		$annotTable->addRow(($$annotNames[$i-1],
				     $$logHash{'annotStats'}[$i]{'readsAligned'},
				     sprintf("%.2f",$$logHash{'annotStats'}[$i]{'cpuTime'})));
	}
	$fileTable->setClass('tableBlue');
    $fileTable->setWidth(600);
    $fileTable->addRow(('Files'));
    $fileTable->addRow('<a href="miR.Counts.csv">miRNA Read Counts</a>');
    $fileTable->addRow('<a href="miR.RPM.csv">miRNA RPM</a>');
    $fileTable->addRow('<a href="mapped.csv">miRNA mapped</a>');
    $fileTable->addRow('<a href="unmapped.csv">miRNA unmapped</a>');
    if($isomirDiff){
        $fileTable->addRow('<a href="isomirs.csv">miRNA isomirs</a>');
        $fileTable->addRow('<a href="isomirs.samples.csv">miRNA isomirs samples</a>');
    }
	open $fh, ">", $filename;
	print $fh htmlHeader();
	print $fh "<h1>miRge Results</h1>\n<h2>Sample Result(s)</h2>\n";
	print $fh $quantTable, "\n<br>\n<h2>Downloadable Files</h2>\n";
    print $fh $fileTable, "\n<br>\n<h2>Annotation summary of unique sequences from sample set and processing time</h2>\n";
	print $fh $annotTable, "\n</body>\n</html>\n";
	close $fh;
}


sub writeDataToCSV {
	my $mappedFile = $outputPath.'/mapped.csv';
	my $isomirFile = $outputPath.'/isomirs.csv';
	my $isomirSampleFile = $outputPath.'/isomirs.samples.csv';
	my $unmappedFile = $outputPath.'/unmapped.csv';
	my $mirRPMFile = $outputPath.'/miR.RPM.csv';
	my $mirCountsFile = $outputPath.'/miR.Counts.csv';
	my $seqKey;
	my $mirKey;
	my $fh;
	my $i;

	open $fh, ">", $mappedFile;
	print $fh "uniqueSequence, annotFlag, $$annotNames[0], $$annotNames[1], $$annotNames[2], $$annotNames[3], $$annotNames[4]";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
		print $fh ", $sampleNames[$i]";
	}
	print $fh "\n";
	# a hash to compare isomirs to their parent mirnas across samples
	# The hash is structured as such:
	# {mirna: {"isomir": {sequence1: [sample 1, .., sample n], sequence2: ...},
	#	   "mirna":  {sequence1: [sample 1, .., sample n], sequence2: ... }}}
	# 
	my %isomirHash = {};
	foreach $seqKey (keys %{$seqHash}) {
		my @entry = ($seqKey);
		my @isomirVals = ($seqKey);
		if ($$seqHash{$seqKey}{'annot'}[0]>0) {
			my $isomir = $$seqHash{$seqKey}{'annot'}[5];
			my $mirna = $$seqHash{$seqKey}{'annot'}[1];
			my $key;
			my $key2;
			if($isomir){
				$key = $isomir;
				$key2 = 'isomirs';
			}
			else{
				$key = $mirna;
				$key2 = 'mirnas';
			}
			# get rid of the SNP annotation for isoMIRs
			if ($key =~ /(.+?)\.SNP/){
				$key = $1;
			}
			unless(exists $isomirHash{$key}){
				$isomirHash{$key} = {};
				$isomirHash{$key}{'mirnas'} = ();
				$isomirHash{$key}{'isomirs'} = {};
			}
			$isomirHash{$key}{$key2}{$seqKey} = ();
			for ($i=0;$i<6;$i++) {
				push(@entry, $$seqHash{$seqKey}{'annot'}[$i]);
			}
			for ($i=0;$i<scalar(@sampleFiles);$i++) {
				my $readCount = $$seqHash{$seqKey}{'quant'}[$i] // 0;
				push(@entry, $readCount);
				push(@{$isomirHash{$key}{$key2}{$seqKey}}, $readCount);
			}
			push(@entry, "\n");
			print $fh join(', ', @entry);
		}
	}
	close $fh;

	if($isomirDiff){
		open $fh, ">", $isomirFile;
		my $sh;
		open $sh, ">", $isomirSampleFile;
		print $sh "miRNA";
		print $fh "miRNA, sequence";
		for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleFiles[$i]";
			print $sh ", $sampleFiles[$i] isomir+miRNA Entropy";
			print $sh ", $sampleFiles[$i] % Canonical Sequence";
			print $sh ", $sampleFiles[$i] Canonical RPM";
			print $sh ", $sampleFiles[$i] Top Isomir RPM";
		}
		print $fh ", Entropy";
		print $fh "\n";
		print $sh "\n";
		my $miRNA;
		my $miRNASeq;
		foreach $miRNA (keys %{isomirHash}) {
			my %sampleIsomirs;
			my $samplemiRNAs = ();
			for ($i=0;$i<scalar(@sampleFiles);$i++) {
				$sampleIsomirs{$i} = ();
				@{$samplemiRNAs}[$i] = 0;
			}
			foreach $miRNASeq (keys %{$isomirHash{$miRNA}{'mirnas'}}) {
				my @sampleArray = @{$isomirHash{$miRNA}{'mirnas'}{$miRNASeq}};
				for ($i=0;$i<scalar(@sampleArray);$i++) {
					@{$samplemiRNAs}[$i] += $sampleArray[$i];
				}
			}
			foreach $miRNASeq (keys %{$isomirHash{$miRNA}{'isomirs'}}) {
				my @entry = ($miRNA);
				push(@entry, $miRNASeq);
				my @sampleArray = @{$isomirHash{$miRNA}{'isomirs'}{$miRNASeq}};
				for ($i=0;$i<scalar(@sampleArray);$i++) {
					push(@{$sampleIsomirs{$i}}, $sampleArray[$i]);
				}
				my $entropy = calcEntropy(\@sampleArray);
				my $maxEntropy = log(scalar(@sampleArray))/log(2);
				if($maxEntropy == 0){
					$entropy = "NA";
				}
				else{
					$entropy = $entropy/$maxEntropy;
				}
				for ($i=0;$i<scalar(@sampleArray);$i++) {
					$sampleArray[$i] *= 1000000.0/$$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'};
				}
				push(@entry, @sampleArray);
				push(@entry, $entropy);
				push(@entry, "\n");
				print $fh join(',', @entry);
			}
			my @isomirOut = ($miRNA);
			for (my $sampleLane=0;$sampleLane<scalar(@sampleFiles);$sampleLane++) {
				my $rpmFactor = 1000000.0/$$logHash{'quantStats'}[$sampleLane]{'mirnaReadsFiltered'};
				my $sampleEntropy = calcEntropy(\@{$sampleIsomirs{$sampleLane}});
				my $topIsomir = max(@{$sampleIsomirs{$sampleLane}})*$rpmFactor;
				my $isomirSum = sumArray(\@{$sampleIsomirs{$sampleLane}})*$rpmFactor;
				push(@{$sampleIsomirs{$sampleLane}}, @{$samplemiRNAs}[$sampleLane]);
				my $sampleEntropyWithmiRNA = calcEntropy(\@{$sampleIsomirs{$sampleLane}});
				my $miRNARPM = $$samplemiRNAs[$sampleLane]*$rpmFactor;

	                        my $maxEntropy = scalar(@{$sampleIsomirs{$sampleLane}});
         	                if($maxEntropy>1){
                 	               $sampleEntropyWithmiRNA = $sampleEntropyWithmiRNA/(log($maxEntropy)/log(2));
	                        }
        	                else{
                	               $sampleEntropyWithmiRNA = "NA";
	                        }

				push(@isomirOut, $sampleEntropyWithmiRNA);

				my $combined = $miRNARPM + $isomirSum;
				if($combined > 0){
					push(@isomirOut, 100*$miRNARPM / $combined);
				}
				else{
					push(@isomirOut, "NA");
				}
				push(@isomirOut, $miRNARPM);
				push(@isomirOut, $topIsomir);
			}
			push(@isomirOut, "\n");
			print $sh join(',', @isomirOut);
		}
		close $fh;
		close $sh;
	}
	
	open $fh, ">", $unmappedFile;
	print $fh "uniqueSequence, annotFlag, $$annotNames[0], $$annotNames[1], $$annotNames[2], $$annotNames[3], $$annotNames[4]";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleNames[$i]";
		}
	print $fh "\n";
	foreach $seqKey (keys %{$seqHash}) {
		if ($$seqHash{$seqKey}{'annot'}[0]==0) {
			print $fh "$seqKey";
			for ($i=0;$i<6;$i++) {
				print $fh ", $$seqHash{$seqKey}{'annot'}[$i]";
			}
			for ($i=0;$i<scalar(@sampleFiles);$i++) {
				print $fh ", ", $$seqHash{$seqKey}{'quant'}[$i] // 0;
			}
			print $fh "\n";
		}
	}
	close $fh;
	
	open $fh, ">", $mirCountsFile;
	print $fh "miRNA";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleNames[$i]";
		}
	print $fh "\n";
	
	print $fh "miRNAtotal";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'}";
		}
	print $fh "\n";

	foreach $mirKey (sort(keys %{$mirHash})) {
		print $fh "$mirKey";
		for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", ", $$mirHash{$mirKey}{'quant'}[$i] // 0;
		}
		print $fh "\n";
	}
	close $fh;
	
	open $fh, ">", $mirRPMFile;
	print $fh "miRNA";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleNames[$i]";
		}
	print $fh "\n";

	foreach $mirKey (sort(keys %{$mirHash})) {
		print $fh "$mirKey";
		for ($i=0;$i<scalar(@sampleFiles);$i++) {
			if($$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'}){
				print $fh ", ", 1000000*$$mirHash{$mirKey}{'quant'}[$i]/$$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'};
			}
			else{
				print $fh ", ", 0;
			}
		}
		print $fh "\n";
	}
	close $fh;
}

sub htmlHeader() {
	return q(<html>
<head>
<style>
.tableBlue {
	margin:0px;padding:0px;
	box-shadow: 5px 5px 5px #888888;
	border:1px solid #000000;
	
	-moz-border-radius-bottomleft:0px;
	-webkit-border-bottom-left-radius:0px;
	border-bottom-left-radius:0px;
	
	-moz-border-radius-bottomright:0px;
	-webkit-border-bottom-right-radius:0px;
	border-bottom-right-radius:0px;
	
	-moz-border-radius-topright:0px;
	-webkit-border-top-right-radius:0px;
	border-top-right-radius:0px;
	
	-moz-border-radius-topleft:0px;
	-webkit-border-top-left-radius:0px;
	border-top-left-radius:0px;
}
.tableBlue table{
    border-collapse: collapse;
        border-spacing: 0;
	width:100%;
	height:100%;
	margin:0px;padding:0px;
}
.tableBlue tr:last-child td:last-child {
	-moz-border-radius-bottomright:0px;
	-webkit-border-bottom-right-radius:0px;
	border-bottom-right-radius:0px;
}
.tableBlue table tr:first-child td:first-child {
	-moz-border-radius-topleft:0px;
	-webkit-border-top-left-radius:0px;
	border-top-left-radius:0px;
}
.tableBlue table tr:first-child td:last-child {
	-moz-border-radius-topright:0px;
	-webkit-border-top-right-radius:0px;
	border-top-right-radius:0px;
}
.tableBlue tr:last-child td:first-child{
	-moz-border-radius-bottomleft:0px;
	-webkit-border-bottom-left-radius:0px;
	border-bottom-left-radius:0px;
}
.tableBlue tr:hover td{
	background-color:#ffffff;
}
.tableBlue td{
	vertical-align:middle;
	background-color:#edf6ff;
	border:1px solid #000000;
	border-width:0px 1px 1px 0px;
	text-align:right;
	padding:7px;
	font-size:12px;
	font-family:Arial;
	font-weight:normal;
	color:#000000;
}
.tableBlue tr:last-child td{
	border-width:0px 1px 0px 0px;
}
.tableBlue tr td:last-child{
	border-width:0px 0px 1px 0px;
}
.tableBlue tr:last-child td:last-child{
	border-width:0px 0px 0px 0px;
}
.tableBlue tr:first-child td{
	background:-o-linear-gradient(bottom, #005fbf 5%, #003f7f 100%);
	background:-webkit-gradient( linear, left top, left bottom, color-stop(0.05, #005fbf), color-stop(1, #003f7f) );
	background:-moz-linear-gradient( center top, #005fbf 5%, #003f7f 100% );
	filter:progid:DXImageTransform.Microsoft.gradient(startColorstr="#005fbf", endColorstr="#003f7f");	background: -o-linear-gradient(top,#005fbf,003f7f);
	background-color:#005fbf;
	border:0px solid #000000;
	text-align:center;
	border-width:0px 0px 1px 1px;
	font-size:14px;
	font-family:Arial;
	font-weight:bold;
	color:#ffffff;
}
.tableBlue tr:first-child:hover td{
	background:-o-linear-gradient(bottom, #005fbf 5%, #003f7f 100%);
	background:-webkit-gradient( linear, left top, left bottom, color-stop(0.05, #005fbf), color-stop(1, #003f7f) );
	background:-moz-linear-gradient( center top, #005fbf 5%, #003f7f 100% );
	filter:progid:DXImageTransform.Microsoft.gradient(startColorstr="#005fbf", endColorstr="#003f7f");	background: -o-linear-gradient(top,#005fbf,003f7f);
	background-color:#005fbf;
}
.tableBlue tr:first-child td:first-child{
	border-width:0px 0px 1px 0px;
}
.tableBlue tr:first-child td:last-child{
	border-width:0px 0px 1px 1px;
}

.thumbnail1  {
	position: relative;
}
.thumbnail1:hover {
	text-decoration: none;
}
.thumbnail1 span { /*CSS for enlarged image*/
	position: absolute;
	background-color: #ffffff;
	padding: 5px;
	border: 4px solid black;
	visibility: hidden;
	color: black;
	text-decoration: none;
}
.thumbnail1:hover span { /*CSS for enlarged image on hover*/
	visibility: visible; 
	z-index: 1;
	top: 0px;
	left: -625px;
}

.thumbnail2  {
	position: relative;
}
.thumbnail2:hover {
	text-decoration: none;
}
.thumbnail2 span { /*CSS for enlarged image*/
	position: absolute;
	background-color: #ffffff;
	padding: 5px;
	border: 4px solid black;
	visibility: hidden;
	color: black;
	text-decoration: none;
}
.thumbnail2:hover span { /*CSS for enlarged image on hover*/
	visibility: visible; 
	z-index: 1;
	top: 0px;
	left: -625px;
}

h1 { font: bold 2.0em Arial; color: #003f7f; text-shadow: 0 0 0.1em #a0a0a0, 0 0 0.1em #a0a0a0; text-align: center}
h2 { font: bold 1.5em Arial; color: #005fbf}

</style>
</head>
);
}


__END__

=head1 SYNOPSIS

perl miRge.pl [--help] [--version] [--adapter none|illumina|ion|sequence] [--species human|mouse|custom_name] [--CPU #] --SampleFiles sample1.fastq,sample2.fastq,...

Examples:

	perl miRge.pl --help
	perl miRge.pl --adapter illumina --species human --SampleFiles seq_file.fastq
	perl miRge.pl --adapter ion --species human --CPU 8 --diff-isomirs --SampleFiles human1.fastq,human2.fastq,human3.fastq
	perl miRge.pl --species mouse --diff-isomirs --phred64 --bowtie /usr/local/bin/bowtie --SampleFiles sample1.fastq,sample2.fastq
	

=head1 OPTIONS

miRge.pl takes the following arguments:

=over 4

=item Required Parameters:
    
=over 4

=item --SampleFiles

						Provide a comma-separated list with no intervening spaces of fastq of fastq.gz
						formatted files containing all of the reads for each individual sample.
						Example: sample1.fastq,sample2.fastq,sample3.fastq
						
=item --species human|mouse|custom_name
										
						Specify which reference species should be used. Used to align with 
						miRge provided references. Custom reference files built with miRge-build.pl
						can also be referenced here by their designated species name.
						
=back

=item Optional Parameters:

=over 4

=item --help

						Displays this usage message.
						
=item --version

						Displays the version of miRge being used.

=item --adapter none|illumina|ion|sequence

						Run adapter removal and quality filtering via cutadapt.
						default: none
						illumina: TGGAATTCTCGGGTGCCAAGGAACTCCAG (TruSeq small RNA kit)
						ion: remove first 11 base pairs (as per miRQC protocol)
						sequence: Provide your own adapter sequence, such as --adapter TGGAATTCTC


=item --CPU #
									
						Specify The number of processors to use for trimming, qc, and alignment.
						default: 1
						
=item --diff-isomirs

						Will output two additional files, xxx.isomirs.samples.csv and xxx.isomirs.csv,
						which are the entropy of isomirs as compared with their canonical miRNAs, and the
						ratio of small RNA species which are canonical. These can be used to flag miRNAs
						that are predominately non-canonical isomiRs and may represent sequencing errors
						and can be used to identify isomir patterns that deviate between specimens.
						default: not generated
						
=item --phred64                                 

						Input fastq files are in phred64 format. By default miRge attempts to inferr the
						phred score from the first 1000 reads of each file. If this fails to detect the phred
						score, this can be used to force phred64.

=item --bowtie                                  

						The path to the system's bowtie binary
						

=item --outputDir

						Optional output path to write analysis results
						
=back


=head2 AUTHORs

=item Alex Baras, E<lt>baras@jhmi.eduE<gt>.

=item Marc Halushka, E<lt>mhalush1@jhmi.eduE<gt>.

=item Chris Mitchell E<lt>cmitch48@jhmi.eduE<gt>.

=item Jason Myers, E<lt>jason_myers@urmc.rochester.eduE<gt>.


=head2 COPYRIGHT

GPL v3

=head2 DATE

01-May-2014

=cut
