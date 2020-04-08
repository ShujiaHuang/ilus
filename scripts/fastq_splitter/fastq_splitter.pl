#!/usr/bin/perl -w
#
#  FASTQ Splitter  -  a script for partitioning a FASTQ file into pieces
#
#  Version 0.2.1 (Jul 31, 2018)
#
#  Copyright (c) 2018 Shujia Huang
#
use Getopt::Long;
use strict;
use File::Basename;


my $start_time = time;

my ($opt_n_parts,$outdir,$opt_version,$uncompress,$opt_help);
my $opt_check = 1;
GetOptions("n-parts=i"   => \$opt_n_parts,
           "outdir=s"    => \$outdir,
           "uncompress"  => \$uncompress,
           "check"       => \$opt_check,
           "version"     => \$opt_version,
           "help"        => \$opt_help) or die "Can't parse command line arguments.\n";

sub show_version {
    print q{FASTQ Splitter 0.2.1 Copyright (c) 2018 Shujia Huang}, "\n";
}

sub show_help {
    print q{Usage: fastq-splitter.pl [options] <fastq file> ...
Options:
    --n-parts <INT>      - Divide into <INT> parts
    --outdir <STR>       - Output diretory. 
    --uncompress         - Uncompress output files.
    --check              - Check FASTQ correctness(1=>check, 0=>not check). [1]
    --version            - Show version.
    --help               - Show help.
},"\n";
}

if ($opt_version) { show_version(); }
if ($opt_help) { show_help(); }
if ($opt_help or $opt_version) {
    # show help message
    exit;
}

if (!defined($outdir)) {
    die "Error: Splitting directory is not specified ('--outdir' option)\n\n";
}

if (!defined($opt_n_parts)) {
    if (!$opt_help and !$opt_version) {
        show_version();
        show_help();
    }
    print STDERR "Error: Divide number is not specified ('--n-parts' option)\n\n";
    exit;
}

die "File for splitting is not specified\n" if !@ARGV;

if (defined($opt_n_parts) and $opt_n_parts <= 0) {
    die "Non-positive number of parts\n";
}

my $n_parts = defined($opt_n_parts) ? $opt_n_parts: 0;

if ($n_parts <= 1) {
    print "WARING: Divide into $n_parts parts, means we don't have to divide anything.\n\n";
    exit(0);
}

# globle value record line number for locating error info.
my $line_count;

foreach my $infile (@ARGV) {
    $line_count = 0;
    print "Loading $infile";
    split_file($infile, $n_parts, $outdir);
}
 
my $end_time = time;
my $elapsed_time = $end_time - $start_time;
print "** Fastq splitter is all done, elapsed $elapsed_time second", (($elapsed_time==1)?'':'s'), " elapsed **\n";

############################################################################################
####################################### Sub function #######################################
sub split_file {

    my ($infile, $n_parts, $outdir) = @_;
    if (!-e $infile or !-f $infile) {
        print "[Error] Can't find file: \"$infile\"\n";
        return;
    }

    my ($fqname, $ext);
    if ($infile =~ m/\.fastq.gz$/) {
        $fqname = basename($infile, ".fastq.gz");
        $ext = "fastq";
    } elsif ($infile =~ m/\.fastq$/) {
        $fqname = basename($infile, ".fastq");
        $ext = "fastq";
    } elsif ($infile =~ m/\.fq.gz$/) {
        $fqname = basename($infile, ".fq.gz");
        $ext = "fq";
    } elsif ($infile =~ m/\.fq$/) {
        $fqname = basename($infile, ".fq");
        $ext = "fq";
    } else {
        die "Error: Input files are not fastq files(*.fq.gz/*.fq/*.fastq.gz/*.fastq)\n";
    }
    $ext .= ".gz" if !$uncompress;

    # Get fastq size and read number
    my ($total_read_num, $total_base_num) = get_file_size($infile);
    my $time = time;
    my $elapsed_time = $time - $start_time;
    print ": $total_read_num sequences, $total_base_num bp, elapsed $elapsed_time seconds elapsed\n";
    print " => dividing into $n_parts parts\n";

    # split fastq
    open(my $IN, ($infile =~ m/\.gz/) ? "gzip -dc $infile |": $infile) or die "Error: Can't open file $infile\n";

    my $num_len = length($n_parts);
    my $read_num_per_file = ($total_read_num % $n_parts == 0) ? $total_read_num / $n_parts : int($total_read_num / $n_parts) + 1;
    for (my $part = 1; $part <= $n_parts; ++$part) {

        my $part_file = sprintf("%s.%0*d.%s",$fqname, $num_len, $part, $ext);
        my $out_sub_file = join("/", $outdir, $part_file);

        print "$out_sub_file\n";

        my $written = 0;
        my ($name_line, $read_seq, $plus_line, $read_qual);
        
        # Compress file by pigz is much much much better than gzip
        open my $OUT, ($out_sub_file =~ m/\.gz$/) ? "| pigz > $out_sub_file": ">$out_sub_file" or die "Error: Cannot write to $out_sub_file\n";
        while(($written < $read_num_per_file) and (!eof $IN)){

            # my ($read_len, $read) = get_fastq_read($IN, $opt_check);
            chomp($name_line=<$IN>);
            chomp($read_seq=<$IN>);
            chomp($plus_line=<$IN>);
            chomp($read_qual=<$IN>);

            my $read = join("\n", ($name_line, $read_seq, $plus_line, $read_qual));
            print $OUT "$read\n";
            ++$written;
        }
        close $OUT;

    }
    close $IN;
}

sub get_file_size {

    my ($file) = @_;
    open(my $IN, ($file =~ m/\.gz/) ? "gzip -dc $file | ": $file) or die "Error: Can't open file $file\n";

    my ($name_line, $read_seq, $plus_line, $read_qual);
    my ($total_read_num, $total_base_num) = (0, 0);
    $line_count = 0;
    while (!eof $IN) {

        # my ($read_len, $read) = get_fastq_read($IN, 1);
        # $total_base_num += $read_len;

        chomp($name_line=<$IN>);
        chomp($read_seq=<$IN>);
        chomp($plus_line=<$IN>);
        chomp($read_qual=<$IN>);
        $line_count += 4;

        ++$total_read_num;
        $total_base_num += length($read_seq);
    }
    close $IN;

    return ($total_read_num, $total_base_num);
}

sub get_fastq_read {
    my ($IN, $check) = @_;
    my $truncated = "Error: Incomplete FASTQ entry at the end -> looks like truncated input!\n";

    if (eof $IN) { return (0, ''); }
    chomp(my $name_line=<$IN>);
    ++$line_count;
    if (substr($name_line, 0, 1) ne '@') {
        die "Error parsing line $line_count: FASTQ entry does not start with '\@':\n$name_line";
    }

    if (eof $IN) { die $truncated; }
    chomp(my $read_seq = <$IN>);
    ++$line_count;
    my $read_len = length($read_seq);

    if (eof $IN) { die $truncated; }
    chomp(my $plus_line = <$IN>);
    ++$line_count;
    if (substr($plus_line,0,1) ne '+') {
        die "Error parsing line $line_count: Expecting '+', found '$plus_line'";
    }

    if (eof $IN) { die $truncated; }
    chomp(my $read_qual = <$IN>);
    ++$line_count;
    my $qual_len = length($read_qual);
    if ($qual_len != $read_len){
        if (eof($IN)) {
            die $truncated;
        } elsif ($check) {
            die "Error: Misformatted FASTQ entry in input line $line_count: 
            quality length ($qual_len) differs from sequence length ($read_len):\n$read_seq\n$read_qual\n";
        }
    }

    my $read = join("\n", ($name_line, $read_seq, $plus_line, $read_qual));
    return ($read_len, $read);
}
