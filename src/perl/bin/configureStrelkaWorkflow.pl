#!/usr/bin/env perl

=head1 LICENSE

Strelka Workflow Software
Copyright (c) 2009-2013 Illumina, Inc.

This software is provided under the terms and conditions of the
Illumina Open Source Software License 1.

You should have received a copy of the Illumina Open Source
Software License 1 along with this program. If not, see
<https://github.com/downloads/sequencing/licenses/>.

=head1 SYNOPSIS

configureStrelkaWorkflow.pl --tumor FILE --normal FILE --ref FILE --config FILE [options]

This script configures the strelka workflow for somatic variant
calling on matched tumor-normal BAM files. The configuration process
will produce an analysis makefile and directory structure. The
makefile can be used to run the analysis on a workstation or compute
cluster via make/qmake or other makefile compatible process.

=head1 ARGUMENTS

=over 4

=item --tumor FILE

Path to tumor sample BAM file (required)

=item --normal FILE

Path to normal sample BAM file (required)

=item --ref FILE

Path to reference genome fasta (required)

=item --config FILE

Strelka configuration file. Default config files can be found in
${STRELKA_INSTALL_ROOT}/etc/ for both ELAND and BWA alignments. (required)

=back

=head1 OPTIONS

=over 4

=item --output-dir DIRECTORY

Root of the analysis directory. This script will place all
configuration files in the analysis directory, and after configuration
all results and intermediate files will be written within the analysis
directory during a run. This directory must not already
exist. (default: ./strelkaAnalysis)

=back

=cut

use warnings FATAL => 'all';
use strict;

use Carp;
$SIG{__DIE__} = \&Carp::confess;

use Cwd qw(getcwd);
use File::Spec;
use Getopt::Long;
use Pod::Usage;

my $baseDir;
my $libDir;
BEGIN {
    my $thisDir=(File::Spec->splitpath($0))[1];
    $baseDir=File::Spec->catdir($thisDir,File::Spec->updir());
    $libDir=File::Spec->catdir($baseDir,'lib');
}
use lib $libDir;
use Utils;

if(getAbsPath($baseDir)) {
    errorX("Can't resolve path for strelka_workflow install directory: '$baseDir'");
}
my $libexecDir=File::Spec->catdir($baseDir,'libexec');
my $optDir=File::Spec->catdir($baseDir,'opt');


my $scriptName=(File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $cmdline = join(' ',$0,@ARGV);


sub usage() { pod2usage(-verbose => 1,
                        -exitval => 2); }

#
# user configuration:
#

my ($tumorBam, $normalBam, $refFile, $configFile, $outDir);
my $help;

GetOptions( "tumor=s" => \$tumorBam,
            "normal=s" => \$normalBam,
            "ref=s" => \$refFile,
            "config=s" => \$configFile,
            "output-dir=s" => \$outDir,
            "help|h" => \$help) || usage();

usage() if($help);
usage() unless($argCount);


#
# Validate input conditions:
#

sub checkFileArg($$) {
   my ($file,$label) = @_;

   errorX("Must specify $label file") unless(defined($file));
   checkFile($file,$label);
}

checkFileArg($tumorBam,"tumor BAM");
checkFileArg($normalBam,"normal BAM");
checkFileArg($refFile,"reference fasta");
checkFileArg($configFile,"configuration ini");

sub makeAbsoluteFilePaths(\$) {
    my ($filePathRef) = @_;

    my ($v,$fileDir,$fileName) = File::Spec->splitpath($$filePathRef);
    if(getAbsPath($fileDir)) {
        errorX("Can't resolve directory path for '$fileDir' from input file argument: '$$filePathRef'");
    }
    $$filePathRef = File::Spec->catfile($fileDir,$fileName);
}

makeAbsoluteFilePaths($tumorBam);
makeAbsoluteFilePaths($normalBam);
makeAbsoluteFilePaths($refFile);
makeAbsoluteFilePaths($configFile);

# also check for BAM index files:
sub checkBamIndex($) {
    my ($file) = @_;
    my $ifile = $file . ".bai";
    if(! -f $ifile) {
        errorX("Can't find index for BAM file '$file'");
    }
}

checkBamIndex($tumorBam);
checkBamIndex($normalBam);


sub checkFaIndex($) {
    my ($file) = @_;
    my $ifile = $file . ".fai";
    if(! -f $ifile) {
        errorX("Can't find index for fasta file '$file'");
    }
    # check that fai file isn't improperly formatted (a la the GATK bundle NCBI 37 fai files)
    open(my $FH,"< $ifile") || errorX("Can't open fai file '$ifile'");
    my $lineno=1;
    while(<$FH>) {
          chomp;
          my @F=split();
          if(scalar(@F) != 5) {
              errorX("Unexpected format for line number '$lineno' of fasta index file: '$ifile'\n\tRe-running fasta indexing may fix the issue. To do so, run: \"samtools faidx $file\"");
          }
          $lineno++;
    }
    close($FH);
}

checkFaIndex($refFile);


if(defined($outDir)) {
    if(getAbsPath($outDir)) {
        errorX("Can't resolve path for ouput directory: '$outDir'");
    }
} else {
    $outDir=File::Spec->catdir(Cwd::getcwd(),'strelkaAnalysis');
}

if(-e $outDir) {
    errorX("Output path already exists: '$outDir'");
}

if(getAbsPath($baseDir)) {
    errorX("Can't resolve path for strelka install directory: '$baseDir'");
}

my $samtoolsDir = File::Spec->catdir($optDir,'samtools');
checkDir($libexecDir,"strelka libexec");
checkDir($samtoolsDir,"samtools");

my $callScriptName = "callSomaticVariants.pl";
my $filterScriptName = "filterSomaticVariants.pl";
my $finishScriptName = "consolidateResults.pl";
my $callScript = File::Spec->catfile($libexecDir,$callScriptName);
my $filterScript = File::Spec->catfile($libexecDir,$filterScriptName);
my $finishScript = File::Spec->catfile($libexecDir,$finishScriptName);
my $countFasta = File::Spec->catfile($libexecDir,"countFastaBases");
my $samtoolsBin = File::Spec->catfile($samtoolsDir,"samtools");
checkFile($callScript,"somatic variant call script");
checkFile($filterScript,"somatic variant filter script");
checkFile($finishScript,"result consolidation script");
checkFile($countFasta,"fasta scanner");
checkFile($samtoolsBin,"samtools");

#
# Configure bin runs:
#
checkMakeDir($outDir);

#
# Configure bin runs: open and validate config ini
#
my $config = parseConfigIni($configFile);

sub checkConfigKeys($) {
    my ($keyref) = @_;
    for (@$keyref) {
        errorX("Undefined configuration option: '$_'") unless(defined($config->{user}{$_}));
    }
}

# these are the keys we need at configuration time:
my @config_keys = qw(binSize);

# these are additional keys we will need to run the workflow:
# (note we don't check for (maxInputDepth,minTier2Mapq) for back compatibility with older config files)
my @workflow_keys = qw(
minTier1Mapq isWriteRealignedBam ssnvPrior sindelPrior ssnvNoise sindelNoise ssnvNoiseStrandBiasFrac
ssnvQuality_LowerBound sindelQuality_LowerBound isSkipDepthFilters depthFilterMultiple
snvMaxFilteredBasecallFrac snvMaxSpanningDeletionFrac indelMaxRefRepeat
indelMaxWindowFilteredBasecallFrac indelMaxIntHpolLength);

checkConfigKeys(\@config_keys);
checkConfigKeys(\@workflow_keys);


my $binSize = int($config->{user}{binSize});

$config->{derived}{configurationCmdline} = $cmdline;
$config->{derived}{normalBam} = $normalBam;
$config->{derived}{tumorBam} = $tumorBam;
$config->{derived}{refFile} = $refFile;
$config->{derived}{outDir} = $outDir;

#
# Configure bin runs: check for consistent chrom info between BAMs and reference
#
sub getBamChromInfo($) {
    my $file = shift;
    my $cmd = "$samtoolsBin view -H $file |";
    open(my $FH,$cmd) || errorX("Can't open process $cmd");

    my %info;
    my $n=0;
    while(<$FH>) {
        next unless(/^\@SQ/);
        chomp;
        my @F = split(/\t/);
        scalar(@F) >= 3 || errorX("Unexpected bam header for file '$file'");

        my %h = ();
        foreach (@F) {
            my @vals = split(':');
            $h{$vals[0]} = $vals[1];
        }
        $F[1] = $h{'SN'};
        $F[2] = $h{'LN'};

        my $size = int($F[2]);
        ($size > 0) || errorX("Unexpected chromosome size '$size' in bam header for file '$file'");
        $info{$F[1]}{size} = $size;
        $info{$F[1]}{order} = $n;
        $n++;
    }
    close($FH) || errorX("Can't close process $cmd");
    return %info;
}


my %chromInfo = getBamChromInfo($normalBam);
my @chroms = sort { $chromInfo{$a}{order} <=> $chromInfo{$b}{order} } (keys(%chromInfo));

{
    #consistency check:
    my %tumorChromInfo = getBamChromInfo($tumorBam);
    for my $chrom (@chroms) {
        my $ln = $chromInfo{$chrom}{size};
        my $tln = $tumorChromInfo{$chrom}{size};
        my $order = $chromInfo{$chrom}{order};
        my $torder = $tumorChromInfo{$chrom}{order};
        unless(defined($tln) && ($tln==$ln) && ($torder==$order)) {
            errorX("Tumor and normal BAM file headers disagree on chromosome: '$chrom'");
        }
        delete $tumorChromInfo{$chrom};
    }
    for my $chrom (keys(%tumorChromInfo)) {
        errorX("Tumor and normal BAM file headers disagree on chromosome: '$chrom'");
    }
}


my %refChromInfo;
logX("Scanning reference genome");
{
    my $knownGenomeSize=0;
    my $cmd="$countFasta $refFile |";
    open(my $FFH,$cmd) || errorX("Failed to open process '$cmd'");

    while(<$FFH>) {
        chomp;
        my @F = split(/\t/);
        scalar(@F) == 4 || errorX("Unexpected value from '$cmd'");
        $knownGenomeSize += int($F[2]);
        $refChromInfo{$F[1]}{knownSize} = int($F[2]);
        $refChromInfo{$F[1]}{size} = int($F[3]);
    }
    close($FFH) || errorX("Failed to close process '$cmd'");

    #consistency check:
    for my $chrom (@chroms) {
        my $ln = $chromInfo{$chrom}{size};
        my $rln = $refChromInfo{$chrom}{size};
        unless(defined($rln) && ($rln==$ln)) {
            errorX("BAM headers and reference fasta disagree on chromosome: '$chrom'");
        }
        $config->{derived}{"chrom_${chrom}_size"} = $rln;
        $config->{derived}{"chrom_${chrom}_knownSize"} = $refChromInfo{$chrom}{knownSize};
    }
    $config->{derived}{chromOrder} = join("\t",@chroms);

    $config->{derived}{knownGenomeSize} = $knownGenomeSize;
}
logX("Scanning reference genome complete");



#
# Configure bin runs: create directory structure
#
my $resultsDir = File::Spec->catdir($outDir,'results');
checkMakeDir($resultsDir);
if($config->{user}{isWriteRealignedBam}) {
    my $bamDir = File::Spec->catdir($outDir,'realigned');
    checkMakeDir($bamDir);
}
my $chromRootDir = File::Spec->catdir($outDir,'chromosomes');
checkMakeDir($chromRootDir);
for my $chrom (@chroms) {
    my $chromDir = File::Spec->catdir($chromRootDir,$chrom);
    checkMakeDir($chromDir);

    my $chromRef = $chromInfo{$chrom};
    $chromRef->{dir} = $chromDir;
    $chromRef->{binList} = getBinList($chromRef->{size},$binSize);

    my $binRootDir = File::Spec->catdir($chromDir,'bins');
    checkMakeDir($binRootDir);

    for my $binId ( @{$chromRef->{binList}} ) {
        my $binDir = File::Spec->catdir($binRootDir,$binId);
        checkMakeDir($binDir);
    }
}



#
# write run config file:
#
my $runConfigFile;
{
    my $cstr = <<END;
;
; Strelka workflow configuration file
;
; This is an automatically generated file, you probably don't want to edit it. If starting a new run,
; input configuration templates (with comments) can be found in the Strelka etc/ directory.
;
END

    $cstr .= writeConfigIni($config);

    my $configDir = File::Spec->catdir($outDir,'config');
    checkMakeDir($configDir);
    $runConfigFile = File::Spec->catdir($configDir,'run.config.ini');
    open(my $FH,"> $runConfigFile") || errorX("Can't open file '$runConfigFile'");
    print $FH $cstr;
    close($FH);
}



#
# create makefile
#
my $makeFile = File::Spec->catfile($outDir,"Makefile");
open(my $MAKEFH, "> $makeFile") || errorX("Can't open file: '$makeFile'");

my $completeFile = "task.complete";

print $MAKEFH <<ENDE;
# This makefile was automatically generated by $scriptName
#
# Please do not edit.

script_dir := $libexecDir
call_script := \$(script_dir)/$callScriptName
filter_script := \$(script_dir)/$filterScriptName
finish_script := \$(script_dir)/$finishScriptName

config_file := $runConfigFile

analysis_dir := $outDir
results_dir := \$(analysis_dir)/results

ENDE

print $MAKEFH <<'ENDE';

complete_tag := task.complete

finish_task := $(analysis_dir)/$(complete_tag)

get_chrom_dir = $(analysis_dir)/chromosomes/$1
get_chrom_task = $(call get_chrom_dir,$1)/$(complete_tag)
get_bin_task = $(call get_chrom_dir,$1)/bins/$2/$(complete_tag)



all: $(finish_task)
	@$(print_success)


define print_success
echo;\
echo Analysis complete. Final somatic calls can be found in $(results_dir);\
echo
endef


# top level results target:
#
$(finish_task):
	$(finish_script) --config=$(config_file) && touch $@


# chromosome targets:
#
ENDE

for my $chrom (@chroms) {

    print $MAKEFH <<ENDE;
chrom_${chrom}_task := \$(call get_chrom_task,$chrom)
\$(finish_task): \$(chrom_${chrom}_task)
\$(chrom_${chrom}_task):
	\$(filter_script) --config=\$(config_file) --chrom=$chrom && touch \$@

ENDE

}

print $MAKEFH <<ENDE;

# chromosome bin targets:
#
ENDE

for my $chrom (@chroms) {
    for my $bin (@{$chromInfo{$chrom}{binList}}) {

print $MAKEFH <<ENDE;
chrom_${chrom}_bin_${bin}_task := \$(call get_bin_task,$chrom,$bin)
\$(chrom_${chrom}_task): \$(chrom_${chrom}_bin_${bin}_task)
\$(chrom_${chrom}_bin_${bin}_task):
	\$(call_script) --config=\$(config_file) --chrom=$chrom --bin=$bin && touch \$@

ENDE

    }
}


# If the eval function is available, this is the way we could finish
# the makefile without being so verbose but it doesn't look like qmake
# understands this function.

=cut

print $MAKEFH <<ENDE;

chroms := @chroms

ENDE

for my $chrom (@chroms) {
    print $MAKEFH "${chrom}_bins := " . join(" ",@{$chromInfo{$chrom}{binList}}) . "\n";
}

print $MAKEFH <<'ENDE';

define chrom_task_template
chrom_$1_task := $(call get_chrom_task,$1)
$(finish_task): $$(chrom_$1_task)
$$(chrom_$1_task):
	$$(filter_script) --config=$$(config_file) --chrom=$1 && touch $$@
endef

$(foreach c,$(chroms),$(eval $(call chrom_task_template,$c)))


# chromosome bin targets:
#
define chrom_bin_task_template
chrom_$1_bin_$2_task := $(call get_bin_task,$1,$2)
$$(chrom_$1_task): $$(chrom_$1_bin_$2_task)
$$(chrom_$1_bin_$2_task):
	$$(call_script) --config=$$(config_file) --chrom=$1 --bin=$2 && touch $$@
endef

$(foreach c,$(chroms), \
    $(foreach b,$($c_bins),$(eval $(call chrom_bin_task_template,$c,$b))) \
 )

ENDE

=cut



print <<END;


Successfully configured analysis and created makefile '$makeFile'.

To run the analysis locally using make, run:

make -C $outDir

...or:

cd $outDir
make

END

1;

__END__

