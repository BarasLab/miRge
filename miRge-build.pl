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
use Time::HiRes qw/ time sleep/;
use Cwd 'abs_path';
use File::Basename;
use File::Spec;


## Path to program and seqLib folder
my $miRgePath = abs_path($0);
$miRgePath =~ s/\/[^\/]+\.pl/\//;
my $refPath = $miRgePath."miRge.seqLibs/";

if(not -d $refPath){
	mkdir $refPath;
}

my $settings={};
my $help;

GetOptions($settings,('help' => \$help,'species=s', 'bowtie-build=s', 'mirna=s', 'hairpin=s','other=s', 'est=s'));

my $speciesType = $$settings{species}||"none";
my $bowtieIndexBuilder = $$settings{'bowtie-build'}||"bowtie-build";
my $bwtIndexDir = File::Spec->catdir($refPath,$speciesType);

my $mirnaBWT = File::Spec->catdir($bwtIndexDir,"mirna");
my $hairpinBWT = File::Spec->catdir($bwtIndexDir,"hairpin");
my $contBWT = File::Spec->catdir($bwtIndexDir,"other");
my $estBWT = File::Spec->catdir($bwtIndexDir,"est");


pod2usage( -verbose => 1) if( $help );
pod2usage({ -message => "You must provide a species type.", -verbose => 0}) if( $speciesType eq "none" );
pod2usage({ -message => "You must provide a miRNA fasta file.", -verbose => 0}) if( not $$settings{mirna} );
pod2usage({ -message => "You must provide a miRNA hairpin fasta file.", -verbose => 0}) if( not $$settings{hairpin} );
pod2usage({ -message => "You must provide an 'other' fasta file.", -verbose => 0}) if( not $$settings{other} );
pod2usage({ -message => "You must provide an est fasta file.", -verbose => 0}) if( not $$settings{est} );


my $mirge_start_time = time;
print "\nChecking for bowtie-build and building indices ...\n";
my $t = time;
checkBowtie();
$t = getTimeDelta($t,time);
print "Bowtie indices built under species name $speciesType ($t sec).\n";

############################# sub-routines #################################
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
	$mirnaBWT = checkBowtieIndex($mirnaBWT, $bowtieIndexBuilder, 'mirna');
	$hairpinBWT = checkBowtieIndex($hairpinBWT, $bowtieIndexBuilder, 'hairpin');
	$contBWT = checkBowtieIndex($contBWT, $bowtieIndexBuilder, 'other');
	$estBWT = checkBowtieIndex($estBWT, $bowtieIndexBuilder, 'est');
}

__END__

=head1 SYNOPSIS

perl miRge-build.pl [--help] [--species species_name] [--mirna fasta file] [--hairpin fasta file] [--other fasta file] [--est fasta file]

Examples:

	perl miRge-build.pl --species human_38 --mirna hsa_mirna.fa --hairpin hsa_hairpin.fa --other human_snorna.fa --est human_refseq.fa
	

=head1 OPTIONS

miRge-build takes the following arguments:

=over 4

=item Required Parameters:
    
=over 4
						
=item --species species_name
										
						The designated species name this group of fasta files will be referenced by.
						The species name corresponds to the directory under which these files are stored.
						As such, we recommend appending a version number or an identifier
						to keep these unique -- such as 'human_38'.
						
=item --mirna					

						A fasta file consisting of mature miRNA reference sequences to align against.
						
=item --hairpin					

						A fasta file consisting of hairpin miRNA reference sequences to align against.
						
=item --other				

						A fasta file consisting of the additional reference sequences to align against.
						This file aims to identify additional sources of small RNA reads such as tRNAs or rRNAs.
						
=item --est					

						A fasta file consisting of additional background sources of miRNAs.
						
=back

=item Optional Parameters:

=over 4


=item --bowtie-build                                  

						The path to your bowtie-build binary
						

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
