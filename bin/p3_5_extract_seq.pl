use strict;
use Getopt::Long;

sub usage{
die "
Description: extract target fasta seq

Usage: perl $0 [Options]
Options:
     -f      fasta file
     -l      list file of target seq id. format: old_id\\tnew_id[optional]
     -o      out file

Example: perl $0 -f all.fasta -l target.list -o result.fa
";
}

my ($output,$all,$list);
&GetOptions(
        'o=s' => \ $output,
        'f=s' => \ $all,
        'l=s' => \ $list
);

unless ( -e $output || -e $all || -e $list){
usage;
}


####################### Main ############################
open IN, "$all" or die $!;
my %seq;
my $id;
while(<IN>){
chomp;
if ($_ =~ /\>/){
$id = (split(" "))[0];
$id =~ s/\>//g;
}else{
$_ =~ s/\s//g;
$seq{$id} .= $_;
}
}
close IN;

open LIST, "$list" or die $!;
open OUT, ">$output";
while(<LIST>){
chomp;
my @ids = split("\t");
my $id;
#$ids[0] =~ s/\s+//g;
my $o_id =(split(" ",$ids[0]))[0];
if ($#ids==1){
$id = $ids[1];
}else{
$id = $o_id;
}
(exists $seq{$o_id})?(print OUT ">$id\n$seq{$o_id}\n"):(print "$o_id can not found in $all\n");
}
close LIST;
close OUT;
