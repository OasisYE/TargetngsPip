#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use File::Path;

=head1 Name
    MainPip.pl -- The main pipeline of 23Mofang.

=head1 Description

=head1 Version
    Author : Weijian Ye ( yeweijian@23mofang.com )
    Date : 2015-09-14

=head1 Usage
    perl MainPip.pl [option]
    
    [option]
    
    --fq1     <string>  The dir of fq1 file
    --fq2     <string>  The dir of fq2 file
    --sam     <string>  Sample ID
    --out     <string>  The outdir of result
    --region  <string>  The target region file,bed format
    --fastQC  <dir>     The path of fastQC 
    --bwa     <dir>     The path of BWA [/usr/src/bwa-0.7.12/bwa]
    --ref     <string>  The path of Human genome reference [GRCh37]
    --annovar <dir>     The dir of annovar program [/data1/annovar]
    -dbsnp    <string>  The path of dbsnp DB 
    --gatk    <string>  The path of GATK
    --gatkM   <int>     The variants calling module of GATK [1:UnifiedGenotyper 2:Haplotypecaller default[1]]
    --start   <int>     Start from which step[1-4 default[1]]
         1    Fastq Quality statistics
         2    Alignment
         3    Variants Detection
         4    Annotation
         5    Final Stat
    --stop    <int>     end of which step[1-5 default[5]]
    --skip    <str>     Skip some steps between [start,stop]
                        e.g.: 1,2,3 means the pipline will skip step 1-3. default [0], means nothing skip.
    --help              show the usuage

Example : 
    perl MainPip.pl [-option] -fq1 fq1File -fq2 fq2File -sam SampleID --region bedFile -out Outdir

=cut

my ( $help, $fq1, $fq2, $sam, $region, $fastQC, $bwa, $samtools, $gatk, $gatkM, $dbsnp, $ref, $annovar, $start, $stop, $skip, $outdir );
GetOptions(
    "help:s"     => \$help,
    "fq1:s"      => \$fq1,
    "fq2:s"      => \$fq2,
    "sam:s"      => \$sam,
    "region:s"   => \$region,
    "ref:s"      => \$ref,
    "fastQC:s"   => \$fastQC,
    "bwa:s"      => \$bwa,
    "samtools:s" => \$samtools,
    "gatk:s"     => \$gatk,
    "dbsnp:s"    => \$dbsnp,
    "annovar:s"  => \$annovar,
    "start:s"    => \$start,
    "stop:s"     => \$stop,
    "skip:s"     => \$skip,
    "gatkm:s"    => \$gatkM,
    "out:s"      => \$outdir,
	
);
die `pod2text $0` if ( $help or !$fq1 or !$fq2 or !$sam or !$region or !$outdir );
$fastQC ||= "/data1/Shark1.1/Tools/fastqc";
$ref ||= "/data1/chenD/SeqMule/database/human_g1k_v37.fasta";
$bwa ||= "/data1/Shark1.1/Tools/bwa-0.7.12/bwa";
$samtools ||= "/data1/Shark1.1/Tools/samtools-1.2/samtools";
$gatk ||= "/data1/Shark1.1/Tools/gatk/GenomeAnalysisTK.jar";
$dbsnp ||= "/data1/Shark1.1/DB/dbsnp144.vcf";
$annovar ||= "/data1/Shark1.1/Tools/annovar";


$start ||= 1;
$stop ||= 5;
$skip ||= 0;
$gatkM ||= 1;


#++++++++++++++++++++++ Main Process ++++++++++++++++++++++#
my %skip;
for ( split /,/, $skip ) { $skip{$_} = $_; }
print STDERR "---Program\t$0\tstarts --> ".localtime()."\n";
print "\n** The pipeline will run from step $start to step $stop, and skip the $skip steps between them **.\n";
print "** The output path: $outdir\n\n";

# Step1 : QC statistics using FastQC
if ( $start <= 1 && $stop >= 1 ){
	if ( !exists $skip{1} ){
		`$fastQC $fq1 $fq2 --extract -t 12 -o $outdir && echo \"Step1 QC statistics done\" > $outdir/step1.QCstat.log`;
		die "Please check step1 QC statistics!!!" if ( !-e "$outdir/step1.QCstat.log" );
	}
}
die "Stop at step $stop\n" if ( $stop <= 1 );

# Step2 : BWA-MEM alignment
if ( $start <= 2 && $stop >= 2 ){
	if ( !exists $skip{2} ){
		`$bwa mem -M -T 0 -A 1 -B 4 -O 6 -E 1 -L 5 -U 17 -t 8 -R '\@RG\tID:$sam\tPL:Illumina\tSM:$sam' $ref $fq1 $fq2 2> $outdir/BwaAlignmentRunning.log | $samtools view -Sh - |$samtools view -uS - |$samtools sort -m 200000000 - $outdir/$sam && $samtools index $outdir/$sam.bam && echo \"Step2 BWA alignment done\" > $outdir/step2.BWAalignment.log`;
		print "$bwa mem -M -T 0 -A 1 -B 4 -O 6 -E 1 -L 5 -U 17 -t 8 -R '\@RG\tID:$sam\tPL:Illumina\tSM:$sam' $ref $fq1 $fq2 2> $outdir/BwaAlignmentRunning.log | $samtools view -Sh - |$samtools view -uS - |$samtools sort -m 200000000 - $outdir/$sam && $samtools index $outdir/$sam.ba\n ";
		die "Please check Step2 BWA alignment!!!" if ( !-e "$outdir/step2.BWAalignment.log" );
	}
}
die "Stop at step $stop\n" if ( $stop <= 2 );

# Step3 : GATK best practices
if ( $start <= 3 && $stop >= 3 ){
	if ( !exists $skip{3} ){
		#3.1 Target Realignment
		`java -Xmx102400m -jar $gatk -T RealignerTargetCreator -R $ref -I $outdir/$sam.bam -L $region -o $outdir/$sam.RealignerTargetCreator.intervals && echo \"Step3.1 GATK RealignerTargetCreator done\" > $outdir/step3.1.RealignerTargetCreator.log`;
		die "Please check Step3.1 GATK TargetRealignment!!!" if ( !-e "$outdir/step3.1.RealignerTargetCreator.log" );

		#3.2 INDEL Realignment
		`java -Xmx102400m -jar $gatk -T IndelRealigner -R $ref -I $outdir/$sam.bam -targetIntervals $outdir/$sam.RealignerTargetCreator.intervals -o $outdir/$sam.IndelRealignment.bam && echo \"Step3.2 GATK IndelRealignment done\" > $outdir/step3.2.IndelRealignment.log`;
		die "Please check Step3.2 GATK IndelRealignment!!!" if ( !-e "$outdir/step3.2.IndelRealignment.log" );

		if ( $gatkM == 2 ){
		#3.3 Haplotypecaller
			`java -Xmx102400m -jar $gatk -T HaplotypeCaller -R $ref -I $outdir/$sam.IndelRealignment.bam --dbsnp $dbsnp -o $outdir/$sam.vcf && echo \"Step3.3 GATK Haplotypecaller done\" > $outdir/step3.3.Haplotypecaller.log`;
			die "Please check Step3.3 GATK Haplotypecaller!!!" if ( !-e "$outdir/step3.3.Haplotypecaller.log" );
		}elsif ( $gatkM == 1 ){
		# 3.3 UnifiedGenotyper
			`java -Xmx102400m -jar $gatk -T UnifiedGenotyper -R $ref -I $outdir/$sam.IndelRealignment.bam --dbsnp $dbsnp -o $outdir/$sam.vcf && echo \"Step3.3 GATK UnifiedGenotyper done\" > $outdir/step3.3.UnifiedGenotyper.log`;
			die "Please check Step3.3 GATK UnifiedGenotyper!!!" if ( !-e "$outdir/step3.3.UnifiedGenotyper.log" );
		}else{
			die "Please imput right gatkM 1:UnifiedGenotyper 2:Haplotypecaller!!!\n";
		}
		# 3.4 GATK VariantFiltration(DP<10 & QD<2.0)
			`java -Xmx102400m -jar $gatk -T VariantFiltration -R $ref -v $outdir/$sam.vcf --filterExpression "QD < 2.0 || DP < 10" --filterName "Filter" -o $outdir/tmp.vcf && less $outdir/tmp.vcf  | perl -e 'while (<>){chomp;if (/^#/){print "\$_\n";next;}\@a=split /\\s+/;if ( \$a[5] =~ /Filter/ ){next;}else{print "\$_\n";}}' > $outdir/$sam.Filter.vcf && rm $outdir/tmp.vcf && echo \"Step3.4 GATK Filter done\" > $outdir/step3.4.Filter.log`;
		die "Please check Step3.4 GATK VariantFiltration!!!" if ( !-e "$outdir/step3.4.Filter.log" );
	}
}
die "Stop at step $stop\n" if ( $stop <= 3 );
	

# Step4 : Annovar Pipeline
if ( $start <= 4 && $stop >= 4 ){
	if ( !exists $skip{4} ){
		`perl $annovar/table_annovar.pl $outdir/$sam.Filter.vcf $annovar/humandb/ -buildver hg19 -out $outdir/$sam.Annovar -remove -protocol refGene,snp138 -operation g,f -nastring . -vcfinput && perl $annovar/table_annovar.pl $outdir/$sam.Annovar.avinput $annovar/humandb/ -buildver hg19 -out $outdir/$sam.AnnovarCSV -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -csvout && echo \"Step4 Annovar Pipeline done\" > $outdir/step4.Annovar.log`;
		die "Please check Step4 Annovar Pipeline!!!" if ( !-e "$outdir/step4.Annovar.log" );
	}
}
die "Stop at step $stop\n" if ( $stop <= 4 );

# Step5 : Final Tidy
if ( $start <= 5 && $stop >= 5 ){
	if ( !exists $skip{5} ){
		#5.1 Make one new dir and multiallele file
		`mkdir -p $outdir/Finalresult && cd $outdir/Finalresult && ln -fs $outdir/$sam.Annovar.hg19_multianno.vcf . && ln -fs $outdir/$sam.AnnovarCSV.hg19_multianno.csv . && less $outdir/$sam.Annovar.hg19_multianno.vcf | perl -e 'while (<>){chomp;next if (/^#/);\@a=split /\\s+/;if ( \$a[4] =~ /\,/ ){print "\$_\n";}}' > $outdir/Finalresult/$sam.Multiallele.txt`;
		#5.2 Build tabix
		`cd $outdir/Finalresult && /data1/Shark1.1/Tools/htslib-1.2.1/bgzip -c $sam.Annovar.hg19_multianno.vcf > $sam.Annovar.hg19_multianno.vcf.gz && /data1/Shark1.1/Tools/htslib-1.2.1/tabix -p vcf $sam.Annovar.hg19_multianno.vcf.gz && rm $sam.Annovar.hg19_multianno.vcf`;
		#5.3 Statistics
		`perl /data1/Shark1.1/Bin/QCtarget.pl -i $outdir/$sam.IndelRealignment.bam -r $region -b /data1/Shark1.1/Lib -o $outdir -plot && cd $outdir/Finalresult && ln -fs $outdir/information.xls .`;
	}
}

print STDERR "All DONE !!!\n";
print STDERR "---Program\t$0\tends  --> ".localtime()."\n";

############################################################
#++++++++++++++++++++++ Sub Function ++++++++++++++++++++++#
############################################################

sub CreateDir {

    my @dir = @_;

    for ( @dir ) {
        system ( "mkdir -p $_" ) if ( ! -d $_ );
    }
}

sub CheckRequireFile { # Make sure that all the reqire files exist. Should check the config file first
    
	my @file = @_;
    for ( @file ) {
        die "[ERROR] The file $_ not exists!\n" if ( ! -f $_ );
    }
}

