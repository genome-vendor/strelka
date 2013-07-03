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

consolidateSomaticVariants.pl [options] | --help

=head2 SUMMARY

Aggregate final results from all chromosomes

=cut

use warnings FATAL => 'all';
use strict;

use Carp;
$SIG{__DIE__} = \&Carp::confess;

use File::Spec;
use File::Temp;
use Getopt::Long;
use Pod::Usage;

my $baseDir;
my $libDir;
my $optDir;
my $vcftDir;
BEGIN {
    my $thisDir=(File::Spec->splitpath($0))[1];
    $baseDir=File::Spec->catdir($thisDir,File::Spec->updir());
    $libDir=File::Spec->catdir($baseDir,'lib');
    $optDir=File::Spec->catdir($baseDir,'opt');
    $vcftDir=File::Spec->catdir($optDir,'vcftools','lib','perl5','site_perl');
}
use lib $libDir;
use Utils;
use lib $vcftDir;
use Vcf;

if(getAbsPath($baseDir)) {
    errorX("Can't resolve path for strelka_workflow install directory: '$baseDir'");
}
$optDir=File::Spec->catdir($baseDir,'opt');


my $scriptName=(File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $cmdline=join(' ',$0,@ARGV);


my $configFile;
my $help;

GetOptions( "config=s" => \$configFile,
            "help|h" => \$help) or pod2usage(2);

pod2usage(2) if($help);
pod2usage(2) unless(defined($configFile));

#
# check fixed paths
#
my $samtoolsBin = File::Spec->catfile($optDir,'samtools','samtools');
checkFile($samtoolsBin,"samtools binary");


#
# read config and validate values
#
checkFile($configFile,"configuration ini");
my $config  = parseConfigIni($configFile);


for (qw(outDir chromOrder)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{derived}{$_}));
}
for (qw(isWriteRealignedBam binSize)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{user}{$_}));
}

my $userconfig = $config->{user};

my @chromOrder = split(/\t/,$config->{derived}{chromOrder});
for my $chrom (@chromOrder) {
    my $chromSizeKey = "chrom_" . $chrom . "_size";
    errorX("Undefined configuration option: '$_'") unless(defined($chromSizeKey));
}

my $outDir = $config->{derived}{outDir};
checkDir($outDir,"output");


my $isWriteRealignedBam = $userconfig->{isWriteRealignedBam};

for my $chrom (@chromOrder) {
    my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);
    checkDir($chromDir,"input chromosome");

    next unless($isWriteRealignedBam);
    my $chromSizeKey = "chrom_" . $chrom . "_size";
    my $binList = getBinList($config->{derived}{$chromSizeKey},$userconfig->{binSize});
    for my $binId (@$binList) {
        my $dir = File::Spec->catdir($chromDir,'bins',$binId);
        checkDir($dir,"input bin");
    }
}



# suffix used for large result file intermediates:
my $itag = ".incomplete";


#
# concatenate vcfs:
#
sub concatenateVcfs($) {
    my $fileName = shift;

    my $is_first = 1;

    my $allFileName = "all." . $fileName;
    my $allFile = File::Spec->catfile($outDir,'results',$allFileName . $itag);
    open(my $aFH,'>',"$allFile")
          || errorX("Failed to open file: '$allFile'");

    # loop over all chroms once to create the header, and one more time for all the data:
    my $headervcf;
    for my $chrom (@chromOrder) {
        my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);
        my $iFile = File::Spec->catfile($chromDir,$fileName);
        checkFile($iFile);

        my $depthKey="maxDepth_${chrom}";

        if($is_first) {
            open(my $iFH,'<',"$iFile")
                || errorX("Failed to open file: '$iFile'");
            $headervcf = Vcf->new(fh=>$iFH);
            $headervcf->parse_header();
            $headervcf->remove_header_line(key=>"cmdline");
            $headervcf->add_header_line({key=>"cmdline",value=>$cmdline});
            $headervcf->remove_header_line(key=>"$depthKey");
            close($iFH);
            $is_first=0;
        }

        {
            open(my $iFH,'<',"$iFile")
                || errorX("Failed to open file: '$iFile'");
            my $vcf = Vcf->new(fh=>$iFH);
            $vcf->parse_header();
            for my $line (@{$vcf->get_header_line(key=>"$depthKey")}) {
                # $line seems to be returned as a length 1 array ref to a hash --  ??!?!??!!
                $headervcf->add_header_line($line->[0]);
            }
            $vcf->close();
            close($iFH);
         }
    }
    print $aFH $headervcf->format_header();
    $headervcf->close();

    for my $chrom (@chromOrder) {
        my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);
        my $iFile = File::Spec->catfile($chromDir,$fileName);

        open(my $iFH,'<',"$iFile")
            || errorX("Failed to open file: '$iFile'");

        my $vcf = Vcf->new(fh=>$iFH);
        $vcf->parse_header();
        print $aFH $_ while(<$iFH>);
    }

    close($aFH);

    # make a second set of files with only the passed variants:
    my $passedFileName = "passed." . $fileName;
    my $passedFile = File::Spec->catfile($outDir,'results',$passedFileName . $itag);
    open(my $pFH,'>',"$passedFile")
          || errorX("Failed to open file: '$passedFile'");

    open(my $arFH,'<',"$allFile")
          || errorX("Failed to open file: '$allFile'");

    while(<$arFH>) {
        chomp;
        unless(/^#/) {
            my @F = split(/\t/);
            next if((scalar(@F)>=7) && ($F[6] ne "PASS"));
        }
        print $pFH "$_\n";
    }

    close($arFH);
    close($pFH);

    my $allFileFinished = File::Spec->catfile($outDir,'results',$allFileName);
    checkMove($allFile,$allFileFinished);

    my $passedFileFinished = File::Spec->catfile($outDir,'results',$passedFileName);
    checkMove($passedFile,$passedFileFinished);
}

concatenateVcfs("somatic.snvs.vcf");
concatenateVcfs("somatic.indels.vcf");


my $bamSuffix = ".realigned.bam";

sub consolidateBam($) {
    my $label = shift;

    my $fileName = $label . $bamSuffix;

    my $reDir = File::Spec->catdir($outDir,'realigned');
    checkMakeDir($reDir);

    my @bamList;
    for my $chrom (@chromOrder) {
        my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);

        my $chromSizeKey = "chrom_" . $chrom . "_size";
        my $binList = getBinList($config->{derived}{$chromSizeKey},$userconfig->{binSize});
        for my $binId (@$binList) {
            my $binDir = File::Spec->catdir($chromDir,'bins',$binId);
            my $rbamFile = File::Spec->catfile($binDir,$fileName);
            checkFile($rbamFile,"bin realigned bam file");

            push @bamList,$rbamFile;
        }
    }

    return unless(scalar(@bamList));

    my $headerFH = File::Temp->new();
    my $getHeaderCmd = "bash -c '$samtoolsBin view -H ".$bamList[0]." > $headerFH'";
    executeCmd($getHeaderCmd);

    my $allFile = File::Spec->catfile($reDir,$fileName . $itag);
    my $cmd="$samtoolsBin merge -h $headerFH $allFile ". join(" ",@bamList);
    executeCmd($cmd);

    my $allFileFinished = File::Spec->catfile($reDir,$fileName);
    checkMove($allFile,$allFileFinished);

    my $indexCmd="$samtoolsBin index $allFileFinished";
    executeCmd($indexCmd);

    # for now don't remove all the bin realignments...
    # unlink(@bamList);
}

if($isWriteRealignedBam) {
    consolidateBam("normal");
    consolidateBam("tumor");
}


1;

__END__

