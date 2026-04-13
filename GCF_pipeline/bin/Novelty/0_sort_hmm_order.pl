use strict;

if (@ARGV != 3){
print "[ERROR] Missing input file, please input 1 standard BGC feature file, 1 wrong BGC feature file and 1 output file\n";
print "[Example] perl $0 <right_hmm_order.tsv> <wrong_hmm_order.tsv> <sorted_hmm_order.tsv>\n";
exit;
}

my @right;
open RIGHT, "$ARGV[0]" or die $!;
while(<RIGHT>){
chomp;
if ($. == 1){
@right = split("\t");
}
next;
}
close RIGHT;

my @n_orders;
open WRONG, "$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]";
while(<WRONG>){
chomp;
my @right_elements;
my $line;
if ($. == 1){
my @tmp = split("\t");
foreach my $r_num (0..$#right){
foreach my $w_num (0..$#tmp){
if ($right[$r_num] eq $tmp[$w_num]){
push @n_orders, $w_num;
push @right_elements, $tmp[$w_num];
}
}
}
}else{
my @tmp = split("\t");
foreach my $num (@n_orders){
push @right_elements, $tmp[$num];
}
}
$line = join("\t", @right_elements);
print OUT "$line\n";
}
close WRONG;
close OUT;

