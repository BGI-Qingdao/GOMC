use strict;

if (@ARGV != 3){
print "[ERROR] missing input file\n";
print "[EXAMPLE] perl $0 <GCF info file> <BGC type file> <output file>\n";
exit;
}

my ($gcf_file,$type_file,$out_file) = @ARGV;

my %types;
open TYPE, $type_file or die $!;
while(<TYPE>){
chomp;
next unless /^\d/;
my ($bgc,$type) = split("\t");
$types{$bgc} = $type unless exists $types{$bgc};
}
close TYPE;

open IN, $gcf_file or die $!;
open OUT, ">$out_file";
print OUT "GCF\tBGC\tMAG\tType\n";
while(<IN>){
chomp;
next unless /^\d/;
my @tmp = split("\t");
$tmp[-1] =~ s/.*;(.*)/\1/g;
my $line = join("\t", @tmp,$types{$tmp[1]});
print OUT "$line\n";
}
close IN;
close OUT;
