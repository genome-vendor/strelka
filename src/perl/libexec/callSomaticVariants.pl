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

callSomaticVariants.pl [options] | --help

=head2 SUMMARY

Run the somatic variant caller for snvs and indels on a single
chromosome bin.

=cut

use warnings FATAL => 'all';
use strict;

use Carp;
$SIG{__DIE__} = \&Carp::confess;

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


my ($chrom, $binId, $configFile);
my $help;

GetOptions(
            "chrom=s" => \$chrom,
            "bin=s" => \$binId,
            "config=s" => \$configFile,
            "help|h" => \$help) or pod2usage(2);

pod2usage(2) if($help);
pod2usage(2) unless(defined($chrom));
pod2usage(2) unless(defined($binId));
pod2usage(2) unless(defined($configFile));



#
# check all fixed paths (not based on commandline arguments):
#
checkDir($baseDir);
checkDir($libexecDir);
checkDir($optDir);

my $strelkaBin=File::Spec->catdir($libexecDir,'strelka2');
checkFile($strelkaBin,"strelka binary");
my $samtoolsBin = File::Spec->catfile($optDir,'samtools','samtools');
checkFile($samtoolsBin,"samtools binary");



#
# read config and validate values
#
checkFile($configFile,"configuration ini");
my $config  = parseConfigIni($configFile);

for (qw(knownGenomeSize tumorBam normalBam refFile outDir)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{derived}{$_}));
}

# we skip maxInputDepth,minTier2Mapq here to allow older config files
for (qw(minTier1Mapq isWriteRealignedBam binSize ssnvPrior sindelPrior
        ssnvNoise sindelNoise ssnvNoiseStrandBiasFrac)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{user}{$_}));
}

my $outDir = $config->{derived}{outDir};
my $binDir = File::Spec->catdir($outDir,'chromosomes',$chrom,'bins',$binId);
checkDir($outDir,"output");
checkDir($binDir,"output bin");


my $tumorBam = $config->{derived}{tumorBam};
my $normalBam = $config->{derived}{normalBam};
my $refFile = $config->{derived}{refFile};
checkFile($tumorBam,"tumor BAM");
checkFile($normalBam,"normal BAM");
checkFile($refFile,"reference");


# pull out some config options for convenience:
my $binSize=$config->{user}{binSize};
my $isWriteRealignedBam=$config->{user}{isWriteRealignedBam};
my $knownGenomeSize = $config->{derived}{knownGenomeSize};


my $begin = (int($binId)*$binSize)+1;
my $end = ((int($binId)+1)*$binSize);
#my $end = $begin+100000;  #debug mode

my $useroptions = $config->{user};

# set previous default value if an older config file is being used:
if(! defined($useroptions->{minTier2Mapq})) {
    $useroptions->{minTier2Mapq} = 5;
}


#
# setup the strelka command-line:
#
my $strelka_base_opts= "-clobber" .
" -filter-unanchored" .
" -min-paired-align-score " . $useroptions->{minTier1Mapq} .
" -min-single-align-score 10" .
" -min-qscore 0" .
" -report-range-begin $begin -report-range-end $end" .
" -samtools-reference '$refFile'" .
" -max-window-mismatch 3 20 -print-used-allele-counts" .
" -bam-seq-name '" . $chrom . "'" .
" -genome-size $knownGenomeSize" .
" -max-indel-size 50" .
" -indel-nonsite-match-prob 0.5" .
" --min-contig-open-end-support 35" .
" --somatic-snv-rate " . $useroptions->{ssnvPrior} .
" --shared-site-error-rate " . $useroptions->{ssnvNoise} .
" --shared-site-error-strand-bias-fraction " . $useroptions->{ssnvNoiseStrandBiasFrac} .
" --somatic-indel-rate " . $useroptions->{sindelPrior} .
" --shared-indel-error-rate " . $useroptions->{sindelNoise} .
" --tier2-min-single-align-score " . $useroptions->{minTier2Mapq} .
" --tier2-min-paired-align-score " . $useroptions->{minTier2Mapq} .
" --tier2-single-align-score-rescue-mode" .
" --tier2-mismatch-density-filter-count 10" .
" --tier2-no-filter-unanchored" .
" --tier2-indel-nonsite-match-prob 0.25" .
" --tier2-include-singleton" .
" --tier2-include-anomalous";


my $somSnvFile='somatic.snvs.unfiltered.vcf';
my $somIndelFile='somatic.indels.unfiltered.vcf';

my $cmd =  "$strelkaBin $strelka_base_opts" .
" -bam-file " . $normalBam .
" --tumor-bam-file " . $tumorBam .
" --somatic-snv-file " . File::Spec->catfile($binDir,$somSnvFile) .
" --somatic-indel-file " . File::Spec->catfile($binDir,$somIndelFile) .
" --variant-window-flank-file 50 " . File::Spec->catfile($binDir,$somIndelFile . '.window');


sub ualignFile($) {
    return File::Spec->catfile($binDir,$_[0] . ".unsorted.realigned.bam");
}
sub alignFile($) {
    return File::Spec->catfile($binDir,$_[0] . ".realigned");
}


if(exists($useroptions->{maxInputDepth}) && ($useroptions->{maxInputDepth} > 0)) {
    $cmd .= " --max-input-depth " . $useroptions->{maxInputDepth};
}


if($isWriteRealignedBam) {
    $cmd .= " -realigned-read-file " . ualignFile("normal") .
            " --tumor-realigned-read-file " . ualignFile("tumor");
}

if(defined($useroptions->{extraStrelkaArguments})){
    my $arg=$useroptions->{extraStrelkaArguments};
    if($arg !~ /^\s*$/) {
        $cmd .= " " . $arg;
    }
}

# this file contains site stats that used to be printed on stdout:
#
$cmd .=
" --report-file " . File::Spec->catfile($binDir,'strelka.stats');

$cmd .=
" >| " . File::Spec->catfile($binDir,'strelka.stdout') .
" 2>| " . File::Spec->catfile($binDir,'strelka.stderr');


executeCmd($cmd,0);


if($isWriteRealignedBam) {
    for my $label (qw(normal tumor)) {
        my $ufile = ualignFile($label);
        if( -f $ufile ) {
            my $afile = alignFile($label);
            my $cmd = "$samtoolsBin sort " . $ufile .  " " . $afile;
            executeCmd($cmd,0);
            unlink($ufile);
        } else {
            logX("Can't find unsorted realigned BAM file: '$ufile'");
        }
    }
}


1;

__END__



