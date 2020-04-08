#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename;

die "USAGE:\nfor SE reads: perl $0 <files number> <FQ> <outdir>\nfor PE reads: perl $0 <files number> <FQ1> <FQ2> <outdir>\n" unless (@ARGV==4 or @ARGV==3);

my $n = shift;
my $outdir = pop;
my @fq=@ARGV;

unless(-d "$outdir")
{
    system("mkdir -p $outdir");
}

if ($n>1){
    my $n_reads = `pigz -dc $fq[0] |wc -l |cut -f1` ;
    $n_reads /= 4 ;

    print STDERR "Total $n_reads pairs in $fq[0]\n";
    my $split_lines = (int( $n_reads / $n) + 1) * 4;

    foreach my $file (@fq){
        my $fqname = basename($file,".fq.gz");
        system("pigz -dc $file | split -d -l $split_lines --filter='pigz > \$FILE.fq.gz' - $outdir/$fqname.") == 0 or die "split failed:  $? " ;
    }
}elsif($n==1){
    foreach my $file (@fq){
        my $fqname = basename($file,".fq.gz");
        system("mv $file $outdir/$fqname.00.fq.gz;");
    }
}else{print "SPLIT number should be greater than 1 !!!\n"}