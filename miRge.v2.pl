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


## Path to programs, by default miRge will use its own copies of the public tools bowtie, cutadapt, and ngsutils within the miRge.seqUtils folder
my $miRgePath = abs_path($0);
$miRgePath =~ s/\/[^\/]+\.pl/\//;
my $cutadaptPath = $miRgePath."miRge.seqUtils/cutadapt-1.7.1/";
my $ngsutilsPath = $miRgePath."miRge.seqUtils/ngsutils/";
my $bowtiePath = $miRgePath."miRge.seqUtils/bowtie/";
my $refPath = $miRgePath."miRge.seqLibs/";

my $settings={};
my $help;
my @sampleFiles;

GetOptions($settings,('help' => \$help,'adapter=s','species=s','CPU=s',
	'SampleFiles=s','isomirCutoff=s', 'bowtie=s', 'mirna=s', 'hairpin=s',
	'contaminants=s', 'est=s', 'cutadapt=s', 'ngsutils=s'));

@sampleFiles = split(',', $$settings{'SampleFiles'});
my $filterFlag = $$settings{adapter}||"none";
my $speciesType = $$settings{species}||0;
my $isomirCutoff = $$settings{isomirCutoff}||0.9;
my $numCPU = $$settings{CPU}||1;
my $cutAdaptBinary = $$settings{cutadapt}||File::Spec->catdir($cutadaptPath, "cutadapt_forked.sh");
my $bowtieBinary = $$settings{bowtie}||File::Spec->catdir($bowtiePath, "bowtie");
my $ngsutilsPath = $$settings{ngsutils}||$ngsutilsPath;
my $mirnaBWT;
my $hairpinBWT;
my $contBWT;
my $estBWT;
# this should be made so the databases are prefixed by the same species name to make it simple
if($speciesType eq 'human'){
	$mirnaBWT = $refPath.$speciesType."_mirna_ebwt/hs_mirna";
	$hairpinBWT = $refPath.$speciesType."_hairpin_ebwt/hs_hairpin";
	$contBWT = $refPath.$speciesType."_snoribturna_ebwt/hs_snoribturna";
	$estBWT = $refPath.$speciesType."_est_ebwt/hs_est";
}
elsif($speciesType eq 'mouse'){
	$mirnaBWT = $refPath.$speciesType."_mirna_ebwt/mmu_mirna";
	$hairpinBWT = $refPath.$speciesType."_hairpin_ebwt/mmu_hairpin";
	$contBWT = $refPath.$speciesType."_snoribturna_ebwt/mmu_snoribturna";
	$estBWT = $refPath.$speciesType."_est_ebwt/mmu_est";
}
else{
	$mirnaBWT = $$settings{mirna};
	$hairpinBWT = $$settings{hairpin};
	$contBWT = $$settings{contaminants};
	$estBWT = $$settings{est};
}

my $seqHash = {}; # key is nucleotide sequence, each record is a hash with
		  # 'annot' a boolean array where element [0] is overall matched status and each element [i] the matched status from each respective round of alignment
		  # 'quant' a integer array of the counts from each sample
my $mirHash = {}; # key is annot string from miRNA library records
my $readLengthHash = {}; # key is length
my $logHash = {}; # two primary keys annot and quant
my $tStamp = int(time);
my $t;
my $annotNames = ['exact miRNA', 'hairpin miRNA', 'non miRNA/mRNA RNA', 'mRNA', 'isomiR miRNA'];

pod2usage( -verbose => 1) if( $help );
pod2usage( -verbose => 1) if( @sampleFiles<1 );
pod2usage( -verbose => 1) if( $numCPU =~ m/\D/ );

print "\nChecking for bowtie and indices ...\n";
$t = time;
checkBowtie();
$t = time-$t;
print "Bowtie and indices present ($t sec).\n";

print "\nStarting quantitation pipeline ...\n";
$t = time;
runQuantitationPipeline();
$t = time-$t;
print "All samples completed ($t sec).\n";

print "\nStarting annotation pipeline ...\n";
$t = time;
runAnnotationPipeline();
$t = time-$t;
print "All annotation cycles completed ($t sec).\n";

print "\nSummarizing and tabulating results ...\n";
$t = time;
summarize();
miRNAmerge();
filter();
#generateGraphs();
writeHtmlReport();
writeDataToCSV();
$t = time-$t;
print "Completed ($t sec)\n";

############################# sub-routines #################################
sub checkBowtieIndex {
	my $path = $_[0];
	if(!(-e $path)){
		die "$path is unable to be found, please check the path.";
	}
	my $builder = $_[1];
	my($filename, $dirs, $suffix) = fileparse($path);
	if($suffix == 'fasta' || $suffix == 'fa' || $suffix == 'fna'){
		my $file = File::Spec->catdir($dirs, $filename);
		my @files = glob("$file.*.ebwt");
		if(scalar(@files)){
			return $file;
		}
		system("$builder $path $file");
		if(-e $file){
			die "Error building $file.";
		}
		return $file;
		
	}
	return -e $path;
}

sub checkBowtie {
	my $bowtieIndex;
	if(-e $bowtieBinary){
		my($filename, $dirs, $suffix) = fileparse($bowtieBinary);
		$bowtieIndex = File::Spec->catdir($dirs, "bowtie-build");
		if(!(-e $bowtieIndex)){
			die "Bowtie-build cannot be found at $bowtieIndex";
		}
	}
	else{
		die "Bowtie cannot be found at $bowtieBinary";
	}
	$mirnaBWT = checkBowtieIndex($mirnaBWT, $bowtieIndex);
	$hairpinBWT = checkBowtieIndex($hairpinBWT, $bowtieIndex);
	$contBWT = checkBowtieIndex($contBWT, $bowtieIndex);
	$estBWT = checkBowtieIndex($estBWT, $bowtieIndex);
}

sub runQuantitationPipeline {
	my $samplePrefix;
	my $cleanedReads;
	my $collapsedReads;
	my $i;
	
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
		print "Processing $sampleFiles[$i] ";
		$$logHash{'quantStats'}[$i]{'filename'} = $sampleFiles[$i];
		$samplePrefix = $sampleFiles[$i];
		$samplePrefix =~ s/\.fastq// ;
		$cleanedReads = $samplePrefix.".trim.fastq";
	
		trimRaw($sampleFiles[$i], $cleanedReads, $i);
		print "cpuTime-trim:$$logHash{'quantStats'}[$i]{'cpuTime-trim'}, ";

		quantReads($cleanedReads, $i);
		print "cpuTime-uniq:$$logHash{'quantStats'}[$i]{'cpuTime-uniq'}\n";
		system("rm $cleanedReads");
		
	}
}

sub trimRaw {
	my $infile = $_[0];
	my $outfile = $_[1];
	my $sampleIndex = $_[2];
	my $fileText;
	my $fh;
	
	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-trim'} = time;
	system("$cutAdaptBinary -q 10 -m 16 -a TGGAATTCTCGGGTGCCAAGGAACTCCAG -e 0.12 --discard-untrimmed -o $outfile $infile > /dev/null");
	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-trim'} = time - $$logHash{'quantStats'}[$sampleIndex]{'cpuTime-trim'};
	
	open $fh, "<", "$infile.log";
	$fileText = join('',<$fh>);
	close $fh;
	system("rm $infile.log");
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

	$$logHash{'quantStats'}[$sampleIndex]{'cpuTime-uniq'} = time - $$logHash{'quantStats'}[$sampleIndex]{'cpuTime-uniq'};
}

sub runAnnotationPipeline {
	my $i;
	my $alignmentStatus;
	my $alignmentResult;	   

	my $lengthFilters = [-26,25,0,0,0];
	my $bwtCmdLines = [
		"$bowtieBinary --threads $numCPU $mirnaBWT -n 0 -f SeqToAnnot.fasta 1>SeqToAnnot.sam 2>SeqToAnnot.log",
		"$bowtieBinary --threads $numCPU $hairpinBWT -n 1 -f SeqToAnnot.fasta 1>SeqToAnnot.sam 2>SeqToAnnot.log",
		"$bowtieBinary --threads $numCPU $contBWT -n 1 -f SeqToAnnot.fasta 1>SeqToAnnot.sam 2>SeqToAnnot.log",
		"$bowtieBinary --threads $numCPU $estBWT -n 0 -f SeqToAnnot.fasta 1>SeqToAnnot.sam 2>SeqToAnnot.log",
		"$bowtieBinary --threads $numCPU $mirnaBWT -f -l 15 -5 1 -3 2 -n 2 SeqToAnnot.fasta 1>SeqToAnnot.sam 2>SeqToAnnot.log"
	]; 
	# -- ALIGNMENT 1 -- length < 26, up to 0 mismatch to miRNA
	# -- ALIGNMENT 2 -- length > 25, up to 1 mismatch to hairpin
	# -- ALIGNMENT 3 -- any length, up to 1 mismatch other RNA
	# -- ALIGNMENT 4 -- any length, up to 0 mismatch to EST
	# -- ALIGNMENT 5 -- any length, up to 2 mismatches with special 5 vs 3 prime considerations to miRNA
	
	for ($i=0; $i<5; $i++) {
		writeSeqToAnnot($$lengthFilters[$i]); # write the records with no matches/annotations and within length criteria (0=none,<0=smaller,>0=greater) to a fasta file to be used by the alignment/matching cmd tools	
		
		# run alignment
		print "Starting Annotation-$$annotNames[$i] ";
		$$logHash{'annotStats'}[$i+1]{'cpuTime'} = time;
		$alignmentStatus = system($$bwtCmdLines[$i]);
		$$logHash{'annotStats'}[$i+1]{'cpuTime'} = time - $$logHash{'annotStats'}[$i+1]{'cpuTime'}; # in seconds
		$alignmentStatus = $alignmentStatus >> 8;
		print "cpuTime:$$logHash{'annotStats'}[$i+1]{'cpuTime'}\n";

		# parse alignment completed without error
		if ($alignmentStatus==0) {
			parseBowtieLog('SeqToAnnot.log',$i+1);
			$alignmentResult = parseAlignment('SeqToAnnot.sam');
			updateAnnotHash($alignmentResult,$i+1);		
		}
		else {
			print "Alignment exited with none-zero status.\n";
			exit;
		}
	
		# remove tmp files
		system("rm SeqToAnnot.fasta SeqToAnnot.sam SeqToAnnot.log");
	}
}

sub writeSeqToAnnot {
	my $lengthFilter = $_[0];
	my $filename = 'SeqToAnnot.fasta';
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
	my $mirList = "$bowtieBinary-inspect -n $estBWT";
		
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
	my $mergeFile;
	my $fh;
	my $line;
	my $i;
	my $j;
	
	if ($speciesType) {
		$mergeFile = $speciesType.'_merges.csv';
	}
	else {
		$mergeFile = 'miRNA_merges.csv';
	}
	if(not -d $refPath){
		mkdir $refPath;
	}
	open $fh, "<", $refPath.$mergeFile;
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
	
	system('mkdir miRge.'.$tStamp.'.graphs');
	
	$lengthMax = max (keys %{$readLengthHash});
	for ($i=0; $i<scalar(@sampleFiles); $i++) {
		undef $graphData;
		for ($j=1; $j<=$lengthMax; $j++) {
			$$graphData[$j] = $$readLengthHash{$j}[$i];
		}
		
		$graph = GD::Graph::bars->new(600,300);
		$countMax = max @{$graphData};
		$graph->set(title => "$sampleFiles[$i] ($$logHash{'quantStats'}[$i]{'trimmedReads'} reads)",
			    x_label => 'Read Length',
			    x_label_position => 0.5,
			    x_label_skip => 5,
			    y_label => 'Counts',
			    y_label_position => 0.5,
			    y_max_value => 10*(1+int($countMax/10)),
			    y_tick_number => 10,
			    borderclrs => undef,
			    dclrs => ['blue']);
		open $fh, ">", 'miRge.'.$tStamp.'.graphs/'.$sampleFiles[$i].'.readDistribution.png';
		print $fh $graph->plot([[0..$lengthMax],$graphData])->png;
		close $fh;
		
		undef $graphData;
		$graph = GD::Graph::hbars->new(600,300);
		$graph->set(title => "$sampleFiles[$i] ($$logHash{'quantStats'}[$i]{'trimmedReads'} reads)",
			    y_label => 'Percentage',
			    y_label_position => 0.5,
			    borderclrs => undef,
			    dclrs => ['blue'],
			    show_values => 1);		
		$graphData = [['miRNA','mRNA','other ncRNA','hairpin','unaligned'],
		[$$logHash{'quantStats'}[$i]{'mirnaReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'mrnaReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'ornaReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'hairpinReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'},$$logHash{'quantStats'}[$i]{'remReads'}/$$logHash{'quantStats'}[$i]{'trimmedReads'}]];
		open $fh, ">", 'miRge.'.$tStamp.'.graphs/'.$sampleFiles[$i].'.readAlignments.png';
		print $fh $graph->plot($graphData)->png;
		close $fh;
	}
}

sub writeHtmlReport {
	my $filename = "miRge.$tStamp.report.html";
	my $fh;
	my $i;
	my $annotTable = new HTML::Table(0,0);
	my $quantTable = new HTML::Table(0,0);	
	
	$quantTable->setClass('tableBlue');
	$quantTable->setWidth(1000);
	$quantTable->addRow(('filename','totalReads','trimmedReads<br>(unique)','mirnaReads / filtered<br>(unique)','hairpinReads','ornaReads','mrnaReads','remReads','composition'));
	for ($i=0; $i<scalar(@sampleFiles); $i++) {
		$quantTable->addRow(($$logHash{'quantStats'}[$i]{'filename'},
				     $$logHash{'quantStats'}[$i]{'totalReads'},
				     '<table><tr></tr><tr><td>'.$$logHash{'quantStats'}[$i]{'trimmedReads'}.'<br>('.$$logHash{'quantStats'}[$i]{'trimmedUniq'}.')</td>'.
				     '<td class="thumbnail1"><img src="miRge.'.$tStamp.'.graphs/'.$sampleFiles[$i].'.readDistribution.png" width="100px" height="50px"><span><img src="miRge.'.$tStamp.'.graphs/'.$sampleFiles[$i].'.readDistribution.png"></span></td></tr></table>',
				     $$logHash{'quantStats'}[$i]{'mirnaReads'}.' / '.$$logHash{'quantStats'}[$i]{'mirnaReadsFiltered'}.'<br>('.$$logHash{'quantStats'}[$i]{'mirnaUniqFiltered'}.')',
				     $$logHash{'quantStats'}[$i]{'hairpinReads'},
				     $$logHash{'quantStats'}[$i]{'ornaReads'},
				     $$logHash{'quantStats'}[$i]{'mrnaReads'},
				     $$logHash{'quantStats'}[$i]{'remReads'},
				     '<table><tr></tr><tr><td class="thumbnail2"><img src="miRge.'.$tStamp.'.graphs/'.$sampleFiles[$i].'.readAlignments.png" width="100px" height="50px"><span><img src="miRge.'.$tStamp.'.graphs/'.$sampleFiles[$i].'.readAlignments.png"></span></td></tr></table>'));
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
	
	open $fh, ">", $filename;
	print $fh htmlHeader();
	print $fh "<h1>miRge 1.1</h1>\n<h2>per sample .fastq file results</h2>\n";
	print $fh $quantTable, "\n<br>\n<h2>annotation of unique sequences from sample set</h2>\n";
	print $fh $annotTable, "\n</body>\n</html>\n";
	close $fh;
}

sub sumArray {
	my @array = $_[0];
	my $arraySum = 0;
	my $i;
	for ($i=0;$i<scalar(@array);$i++){
		$arraySum += $array[$i];
	}
	return $arraySum;
}

sub writeDataToCSV {
	my $mappedFile = "miRge.$tStamp.mapped.csv";
	my $isomirFile = "miRge.$tStamp.isomirs.csv";
	my $unmappedFile = "miRge.$tStamp.unmapped.csv";
	my $mirRPMFile = "miRge.$tStamp.miR.RPM.csv";
	my $mirCountsFile = "miRge.$tStamp.miR.Counts.csv";
	my $seqKey;
	my $mirKey;
	my $fh;
	my $i;

	open $fh, ">", $mappedFile;
	print $fh "uniqueSequence, annotFlag, $$annotNames[0], $$annotNames[1], $$annotNames[2], $$annotNames[3], $$annotNames[4]";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleFiles[$i]";
		}
	print $fh "\n";
	# a hash to compare isomirs to their parent mirnas across samples
	my %isomirHash;
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
		$isomirHash{$i}{'isomirs'} = {};
		$isomirHash{$i}{'mirnas'} = {};
	}
	foreach $seqKey (keys %{$seqHash}) {
		my @entry = ($seqKey);
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
			for ($i=0;$i<6;$i++) {
				push(@entry, $$seqHash{$seqKey}{'annot'}[$i]);
			}
			for ($i=0;$i<scalar(@sampleFiles);$i++) {
				my $readCount = $$seqHash{$seqKey}{'quant'}[$i] // 0;
				push(@entry, $readCount);
				if(exists $isomirHash{$i}{$key}{$key2}){
					$isomirHash{$i}{$key}{$key2} += $readCount;
				}
				else{
					$isomirHash{$i}{$key}{$key2} = $readCount;
				}
			}
			push(@entry, "\n");
			print $fh join(', ', @entry);
		}
	}
	close $fh;

	open $fh, ">", $isomirFile;
	print $fh "uniqueSequence, annotFlag, $$annotNames[0], $$annotNames[1], $$annotNames[2], $$annotNames[3], $$annotNames[4]";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleFiles[$i]";
		}
	print $fh "\n";
	foreach $seqKey (keys %{$seqHash}) {
		my $mirna = $$seqHash{$seqKey}{'annot'}[5];
		if(!$mirna){
			$mirna = $$seqHash{$seqKey}{'annot'}[1];
		}
		for ($i=0;$i<scalar(@sampleFiles);$i++) {
			my $mirnaCount = $isomirHash{$i}{$mirna}{'mirnas'};
			my $isomirCount = $isomirHash{$i}{$mirna}{'isomirs'};
			if(!$mirnaCount){
				$mirnaCount = 0;
			}
			if(!$isomirCount){
				$isomirCount = 0;
			}
			print "$isomirCount, $mirnaCount, $mirna \n";
			if($isomirCount+$mirnaCount > 0){
				my $isomirRatio = $isomirCount/($isomirCount+$mirnaCount);
				if($isomirRatio > $isomirCutoff){
					my @entry = ($seqKey);
					for ($i=0;$i<6;$i++) {
						push(@entry, $$seqHash{$seqKey}{'annot'}[$i]);
					}
					for ($i=0;$i<scalar(@sampleFiles);$i++) {
						push(@entry, $$seqHash{$seqKey}{'quant'}[$i] // 0);
					}
					push(@entry, "\n");
					print $fh join(', ', @entry);
					last;
				}
			}
		}
	}
	close $fh;
	
	open $fh, ">", $unmappedFile;
	print $fh "uniqueSequence, annotFlag, $$annotNames[0], $$annotNames[1], $$annotNames[2], $$annotNames[3], $$annotNames[4]";
	for ($i=0;$i<scalar(@sampleFiles);$i++) {
			print $fh ", $sampleFiles[$i]";
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
			print $fh ", $sampleFiles[$i]";
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
			print $fh ", $sampleFiles[$i]";
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
	left: 120px;
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

h1 { font: bold 2.0em Arial; color: #003f7f; text-shadow: 0 0 0.1em #a0a0a0, 0 0 0.1em #a0a0a0}
h2 { font: bold 1.5em Arial; color: #005fbf}

</style>
</head>
);
}


__END__

=head1 SYNOPSIS

perl miRge.v2.pl [--help] [--man] [--cutadapt none|illumina|ion] [--species human|mouse] [--CPU #] --SampleFiles sample1.fastq,sample2.fastq,...

Examples:

	perl miRge.v1.pl --help
	perl miRge.v1.pl --SampleFiles human1.fastq,human2.fastq
	perl miRge.v1.pl --cutadapt illumina --species mouse --CPU 12 --SampleFiles [sample1.fastq,sample2.fastq]
	

=head2 ARGUMENTS

miRge.v1.pl takes the following arguments:

=over 4

=item Optional Parameters:

=over 4

=item --help

						Displays the usage message.

=item --cutadapt none|illumina|ion

						Run adapter removal and quality filtering via cutadapt.
						default: none
						illumina: TGGAATTCTCGGGTGCCAAGGAACTCCAG (TruSeq small RNA kit)
						ion: remove first 11 base pairs (as per miRQC protocol)

=item --species human|mouse
										
						Specify which reference species should be used. Used to align with 
						miRge provided references.
						default: human

=item --CPU #
									
						Sepcify The number of processors to use for trimming, qc, and alignment.
						default: 1

=back

=item Required Parameters:
    
=over 4

=item --SampleFiles

						Provide a comma-seperated list with no intervening space of fastq
						formatted files	containing all of the reads for each individual sample.

=back

=head2 AUTHORs

=item Alex Baras, E<lt>baras@jhmi.eduE<gt>.

=item Marc Halushka, E<lt>mhalush1@jhmi.eduE<gt>.

=item Jason Myers, E<lt>jason_myers@urmc.rochester.eduE<gt>.

=head2 COPYRIGHT

GPL

=head2 DATE

01-May-2014

=cut
