#!/usr/bin/perl -w
#
#  Version 0.1.1 (Aug 29, 2018)
#
#  Copyright (c) 2018 Shujia Huang
#
use strict;

my ($fafile) = @ARGV;

my %facontent;
open I, $fafile or die "Cannot open file : $fafile\n";

my $seq;
while (<I>) {

    next if /^>/;
    chomp;

    my @seq = split(//);
    for my $b (@seq) {
        $facontent{uc($b)}++;
    }
}
close I;


for my $k (keys %facontent){
    print "$k\t$facontent{$k}\n";
}
