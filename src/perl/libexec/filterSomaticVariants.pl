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

filterSomaticVariants.pl [options] | --help

=head2 SUMMARY

Aggregate and filter the variant caller results by chromosome

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
my $vcftDir;
BEGIN {
    my $thisDir=(File::Spec->splitpath($0))[1];
    $baseDir=File::Spec->catdir($thisDir,File::Spec->updir());
    $libDir=File::Spec->catdir($baseDir,'lib');
    my $optDir=File::Spec->catdir($baseDir,'opt');
    $vcftDir=File::Spec->catdir($optDir,'vcftools','lib','perl5','site_perl');
}
use lib $libDir;
use Utils;
use lib $vcftDir;
use Vcf;


my $scriptName=(File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $cmdline = join(' ',$0,@ARGV);


my ($chrom, $configFile);
my $help;

GetOptions( "chrom=s" => \$chrom,
            "config=s" => \$configFile,
            "help|h" => \$help) or pod2usage(2);

pod2usage(2) if($help);
pod2usage(2) unless(defined($chrom));
pod2usage(2) unless(defined($configFile));



#
# read config and validate values
#
checkFile($configFile,"configuration ini");
my $config  = parseConfigIni($configFile);

my $chromSizeKey = "chrom_" . $chrom . "_size";
my $chromKnownSizeKey = "chrom_" . $chrom . "_knownSize";

for (("outDir", $chromSizeKey, $chromKnownSizeKey)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{derived}{$_}));
}
for (qw(binSize ssnvQuality_LowerBound sindelQuality_LowerBound
        isSkipDepthFilters depthFilterMultiple
        snvMaxFilteredBasecallFrac snvMaxSpanningDeletionFrac
        indelMaxRefRepeat indelMaxWindowFilteredBasecallFrac
        indelMaxIntHpolLength)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{user}{$_}));
}

my $outDir = $config->{derived}{outDir};
my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);
checkDir($outDir,"output");
checkDir($chromDir,"output chromosome");

my $userconfig = $config->{user};

my $binList = getBinList($config->{derived}{$chromSizeKey},$userconfig->{binSize});
for my $binId (@$binList) {
    my $dir = File::Spec->catdir($chromDir,'bins',$binId);
    checkDir($dir,"input bin");
}

#
# parameters from user config file:
#

# minimum passed ssnv_nt Q-score:
my $minQSSNT=$userconfig->{ssnvQuality_LowerBound};
# minimum passed sindel_nt Q-score:
my $minQSINT=$userconfig->{sindelQuality_LowerBound};

#
# filtration parameters from user config file:
#

# skip depth filters for targeted resequencing:
my $isUseDepthFilter=(! $userconfig->{isSkipDepthFilters});
# multiple of the normal mean coverage to filter snvs and indels
my $depthFilterMultiple=$userconfig->{depthFilterMultiple};
# max filtered basecall fraction for any sample
my $snvMaxFilteredBasecallFrac=$userconfig->{snvMaxFilteredBasecallFrac};
# max snv spanning-deletion fraction for any sample
my $snvMaxSpanningDeletionFrac=$userconfig->{snvMaxSpanningDeletionFrac};
# max indel reference repeat length
my $indelMaxRefRepeat=$userconfig->{indelMaxRefRepeat};
# max indel window filter fraction
my $indelMaxWindowFilteredBasecallFrac=$userconfig->{indelMaxWindowFilteredBasecallFrac};
# max indel interupted hompolymer length:
my $indelMaxIntHpolLength=$userconfig->{indelMaxIntHpolLength};


# first we want the normal sample mean chromosome depth:
#
my $filterCoverage;
if($isUseDepthFilter) {
    my $totalCoverage = 0;
    for my $binId (@$binList) {
        my $sFile = File::Spec->catfile($chromDir,'bins',$binId,'strelka.stats');
        checkFile($sFile,"strelka bin stats");
        open(my $sFH, '<', $sFile)
          || errorX("Can't open file: '$sFile' $!");

        my $is_found=0;
        while(<$sFH>) {
            next unless(/^NORMAL_NO_REF_N_COVERAGE\s/);
            my $is_bad = 0;
            if(not /sample_size:\s*(\d+)\s+/) { $is_bad=1; }
            my $ss=$1;
            # leave the regex for mean fairly loose to pick up scientific notation, etc..
            if(not /mean:\s*(\d[^\s]*|nan)\s+/) { $is_bad=1; }
            my $mean=$1;
            errorX("Unexpected format in file: '$sFile'") if($is_bad);

            my $coverage = ( $ss==0 ? 0 : int($ss*$mean) );

            $totalCoverage += $coverage;
            $is_found=1;
            last;
        }
        close($sFH);
        errorX("Unexpected format in file: '$sFile'") unless($is_found);
    }

    my $chromKnownSize = $config->{derived}{$chromKnownSizeKey};
    my $normalMeanCoverage = ($totalCoverage/$chromKnownSize);
    $filterCoverage = $normalMeanCoverage*$depthFilterMultiple;
}


# add filter description to vcf header unless it already exists
# return 1 if filter id already exists, client can decide if this is an error
#
sub add_vcf_filter($$$) {
    my ($vcf,$id,$desc) = @_;
    return 1 if(scalar(@{$vcf->get_header_line(key=>'FILTER', ID=>$id)}));
    $vcf->add_header_line({key=>'FILTER', ID=>$id,Description=>$desc});
    return 0;
}


sub check_vcf_for_sample($$$) {
    my ($vcf,$sample,$file) = @_;
    my $is_found=0;
    for ($vcf->get_samples()) {
        $is_found=1 if($_ eq $sample);
    }
    errorX("Failed to find sample '$sample' in vcf file '$file'") unless($is_found);
}



my $depthFiltId="DP";


# Runs all post-call vcf filters:
#
sub filterSnvFileList(\@$$$) {
    my ($ifiles,$depthFilterVal,$acceptFileName,$isUseDepthFilter) = @_;

    my $baseFiltId="BCNoise";
    my $spanFiltId="SpanDel";
    my $qFiltId="QSS_ref";

    open(my $aFH,'>',$acceptFileName)
      or errorX("Failed to open file: '$acceptFileName'");

    my $is_first=1;
    for my $ifile (@$ifiles) {
        open(my $iFH,'<',"$ifile")
          or errorX("Failed to open file: '$ifile'");
        my $vcf = Vcf->new(fh=>$iFH);

        # run some simple header validation on each vcf:
        $vcf->parse_header();
        check_vcf_for_sample($vcf,'NORMAL',$ifile);
        check_vcf_for_sample($vcf,'TUMOR',$ifile);

        if($is_first) {
            # TODO: update vcf meta-data for chromosome-level filtered files
            #
            $vcf->remove_header_line(key=>"cmdline");
            $vcf->add_header_line({key=>"cmdline",value=>$cmdline});
            if($isUseDepthFilter) {
                $vcf->add_header_line({key=>"maxDepth_$chrom",value=>$depthFilterVal});
                add_vcf_filter($vcf,$depthFiltId,"Greater than ${depthFilterMultiple}x chromosomal mean depth in Normal sample");
            }
            add_vcf_filter($vcf,$baseFiltId,"Fraction of basecalls filtered at this site in either sample is at or above $snvMaxFilteredBasecallFrac");
            add_vcf_filter($vcf,$spanFiltId,"Fraction of reads crossing site with spanning deletions in either sample exceeeds $snvMaxSpanningDeletionFrac");
            add_vcf_filter($vcf,$qFiltId,"Normal sample is not homozygous ref or ssnv Q-score < $minQSSNT, ie calls with NT!=ref or QSS_NT < $minQSSNT");
            print $aFH $vcf->format_header();
            $is_first=0;
        }

        while(my $x=$vcf->next_data_hash()) {
            my $norm=$x->{gtypes}->{NORMAL};
            my $tumr=$x->{gtypes}->{TUMOR};

            my %filters;

            # normal depth filter:
            my $normalDP=$norm->{DP};
            if($isUseDepthFilter) {
                $filters{$depthFiltId} = ($normalDP > $depthFilterVal);
            }

            # filtered basecall fraction:
            my $normal_filt=($normalDP>0 ? $norm->{FDP}/$normalDP : 0);

            my $tumorDP=$tumr->{DP};
            my $tumor_filt=($tumorDP>0 ? $tumr->{FDP}/$tumorDP : 0);

            $filters{$baseFiltId}=(($normal_filt >= $snvMaxFilteredBasecallFrac) or
                                   ($tumor_filt >= $snvMaxFilteredBasecallFrac));

            # spanning deletion fraction:
            my $normalSDP=$norm->{SDP};
            my $normalSpanTot=($normalDP+$normalSDP);
            my $normalSpanDelFrac=($normalSpanTot>0 ? ($normalSDP/$normalSpanTot) : 0);

            my $tumorSDP=$tumr->{SDP};
            my $tumorSpanTot=($tumorDP+$tumorSDP);
            my $tumorSpanDelFrac=($tumorSpanTot>0 ? ($tumorSDP/$tumorSpanTot) : 0);

            $filters{$spanFiltId}=(($normalSpanDelFrac > $snvMaxSpanningDeletionFrac) or
                                   ($tumorSpanDelFrac > $snvMaxSpanningDeletionFrac));

            # Q-val filter:
            $filters{$qFiltId}=(($x->{INFO}->{NT} ne "ref") or
                                ($x->{INFO}->{QSS_NT} < $minQSSNT));

            $x->{FILTER} = $vcf->add_filter($x->{FILTER},%filters);

            print $aFH $vcf->format_line($x);
        }

        $vcf->close();
        close($iFH);
    }

    close($aFH);
}



sub updateA2(\@$) {
    my ($a2,$FH) = @_;
    my $line=<$FH>;
    unless(defined($line)) { errorX("Unexpected end of somatic indel window file"); }
    chomp $line;
    @$a2 = split("\t",$line);
    unless(scalar(@$a2)) { errorX("Unexpected format in somatic indel window file"); }
}



sub filterIndelFileList(\@$$$) {
    my ($ifiles,$depthFilterVal,$acceptFileName,$isUseDepthFilter) = @_;

    my $repeatFiltId="Repeat";
    my $iHpolFiltId="iHpol";
    my $baseFiltId="BCNoise";
    my $qFiltId="QSI_ref";

    open(my $aFH,'>',$acceptFileName)
      or errorX("Failed to open file: '$acceptFileName'");

    my $is_first=1;
    for my $ifile (@$ifiles) {
        open(my $iFH,'<',"$ifile")
          or errorX("Failed to open somatic indel file: '$ifile'");

        my $iwfile = $ifile . ".window";
        open(my $iwFH,'<',"$iwfile")
          or errorX("Failed to open somatic indel window file: '$iwfile'");

        my @a2; # hold window file data for one line in case we overstep...

        my $vcf = Vcf->new(fh=>$iFH);

        # run some simple header validation on each vcf:
        $vcf->parse_header();
        check_vcf_for_sample($vcf,'NORMAL',$ifile);
        check_vcf_for_sample($vcf,'TUMOR',$ifile);

        if($is_first) {
            # TODO: update all vcf meta-data for chromosome-level filtered files
            #
            $vcf->remove_header_line(key=>"cmdline");
            $vcf->add_header_line({key=>"cmdline",value=>$cmdline});
            if($isUseDepthFilter) {
                $vcf->add_header_line({key=>"maxDepth_$chrom",value=>$depthFilterVal});
            }
            $vcf->add_header_line({key=>'FORMAT', ID=>'DP50',Number=>1,Type=>'Float',Description=>'Average tier1 read depth within 50 bases'});
            $vcf->add_header_line({key=>'FORMAT', ID=>'FDP50',Number=>1,Type=>'Float',Description=>'Average tier1 number of basecalls filtered from original read depth within 50 bases'});
            $vcf->add_header_line({key=>'FORMAT', ID=>'SUBDP50',Number=>1,Type=>'Float',Description=>'Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases'});
            if($isUseDepthFilter) {
                add_vcf_filter($vcf,$depthFiltId,"Greater than ${depthFilterMultiple}x chromosomal mean depth in Normal sample");
            }
            add_vcf_filter($vcf,$repeatFiltId,"Sequence repeat of more than ${indelMaxRefRepeat}x in the reference sequence");
            add_vcf_filter($vcf,$iHpolFiltId,"Indel overlaps an interupted homopolymer longer than ${indelMaxIntHpolLength}x in the reference sequence");
            add_vcf_filter($vcf,$baseFiltId,"Average fraction of filtered basecalls within 50 bases of the indel exceeds $indelMaxWindowFilteredBasecallFrac");
            add_vcf_filter($vcf,$qFiltId,"Normal sample is not homozygous ref or sindel Q-score < $minQSINT, ie calls with NT!=ref or QSI_NT < $minQSINT");
            print $aFH $vcf->format_header();
            $is_first=0;
        }

        while(my $x=$vcf->next_data_hash()) {
            $vcf->add_format_field($x,'DP50');
            $vcf->add_format_field($x,'FDP50');
            $vcf->add_format_field($x,'SUBDP50');

            my $norm=$x->{gtypes}->{NORMAL};
            my $tumr=$x->{gtypes}->{TUMOR};

            my $chrom=$x->{CHROM};
            my $pos=int($x->{POS});

            # get matching line from window file:
            while((scalar(@a2)<2) or
                  (($a2[0] le $chrom) and (int($a2[1]) < $pos))) {
                updateA2(@a2,$iwFH);
            }
            unless(scalar(@a2) and ($a2[0] eq $chrom) and (int($a2[1]) == $pos))
                { errorX("Can't find somatic indel window position.\nIndel line: " . $vcf->format_line($x) ); }

            # add window data to vcf record:
            #
            $norm->{DP50} = $a2[2]+$a2[3];
            $norm->{FDP50} = $a2[3];
            $norm->{SUBDP50} = $a2[4];
            $tumr->{DP50} = $a2[5]+$a2[6];
            $tumr->{FDP50} = $a2[6];
            $tumr->{SUBDP50} = $a2[7];

            my %filters;

            # normal depth filter:
            my $normalDP=$norm->{DP};
            if($isUseDepthFilter) {
                $filters{$depthFiltId}=($normalDP > $depthFilterVal);
            }

            # ref repeat
            my $refRep=$x->{INFO}->{RC};
            $filters{$repeatFiltId}=(defined($refRep) and
                                     ($refRep > $indelMaxRefRepeat));

            # ihpol
            my $iHpol=$x->{INFO}->{IHP};
            $filters{$iHpolFiltId}=(defined($iHpol) and
                                    ($iHpol > $indelMaxIntHpolLength));

            # base filt:
            my $normWinFrac=( $norm->{DP50} ? $norm->{FDP50}/$norm->{DP50} : 0 );
            my $tumrWinFrac=( $tumr->{DP50} ? $tumr->{FDP50}/$tumr->{DP50} : 0 );
            $filters{$baseFiltId}=( ($normWinFrac >= $indelMaxWindowFilteredBasecallFrac) or
                                    ($tumrWinFrac >= $indelMaxWindowFilteredBasecallFrac) );

            # Q-val filter:
            $filters{$qFiltId}=(($x->{INFO}->{NT} ne "ref") or
                                ($x->{INFO}->{QSI_NT} < $minQSINT));

            $x->{FILTER} = $vcf->add_filter($x->{FILTER},%filters);

            print $aFH $vcf->format_line($x);
        }

        $vcf->close();
        close($iFH);
        close($iwFH);
    }

    close($aFH);
}



my @ssnvFiles;
my @sindelFiles;
for my $binId (@$binList) {
    my $ssnvFile = File::Spec->catfile($chromDir,'bins',$binId,'somatic.snvs.unfiltered.vcf');
    my $sindelFile = File::Spec->catfile($chromDir,'bins',$binId,'somatic.indels.unfiltered.vcf');
    checkFile($ssnvFile,"bin snv file");
    checkFile($sindelFile,"bin indel file");
    push @ssnvFiles,$ssnvFile;
    push @sindelFiles,$sindelFile;
}

my $ssnvOutFile = File::Spec->catfile($chromDir,"somatic.snvs.vcf");
filterSnvFileList(@ssnvFiles,$filterCoverage,$ssnvOutFile,$isUseDepthFilter);

my $sindelOutFile = File::Spec->catfile($chromDir,"somatic.indels.vcf");
filterIndelFileList(@sindelFiles,$filterCoverage,$sindelOutFile,$isUseDepthFilter);


1;
