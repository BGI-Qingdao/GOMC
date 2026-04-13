use strict;

`mkdir -p $ARGV[0]/Circlize_Input` unless -e "$ARGV[0]/Circlize_Input";
my %phylums;
if ($ARGV[1] eq "T"){
	open TOTAL, "$ARGV[0]/temp/GCF_counts_per_phylim_total_sorted.txt" or die $!;
}else{
	open TOTAL, "$ARGV[0]/temp/GCF_counts_per_phylim_total.txt" or die $!;
}
open CHR, ">$ARGV[0]/Circlize_Input/Chr_len.txt";
print CHR "chr\tstart\tend\n";
open INFO, ">$ARGV[0]/Circlize_Input/Chr_info.txt";
print INFO "Phylum\tNumber\n";
while (<TOTAL>){
	chomp;
	next if $. == 1;
	my @tmp = split("\t");
	$phylums{$tmp[0]} = $.-1;
	my $line = join("\t", $phylums{$tmp[0]},0,$tmp[1]);
	print CHR "$line\n";
	print INFO "$tmp[0]\t$phylums{$tmp[0]}\n";
}
close TOTAL;
close CHR;
close INFO;

open UNIQ, "$ARGV[0]/temp/GCF_counts_per_phylim_uniq.txt" or die $!;
my %uniq;
my %count;
while(<UNIQ>){
	chomp;
	next if $. ==1;
	my @tmp = split("\t");
	$uniq{$phylums{$tmp[0]}} = $tmp[1];
}
close UNIQ;

open LINK_OUT, ">$ARGV[0]/Circlize_Input/Links_input.txt" or die $!;
foreach my $phylum (sort keys %uniq){
	my $line = join("\t", $phylum,0,$uniq{$phylum},$phylum,0,$uniq{$phylum});
	print LINK_OUT "$line\n";
}

if ($ARGV[2] eq "T"){
	open LINK, "$ARGV[0]/temp/link_simplify.txt" or die $!;
}else{
	open LINK, "$ARGV[0]/temp/link.txt" or die $!;
}
my %phylums_u;
while (<LINK>){
	chomp;
	my @tmp = split("\t");
    unless (grep /^$tmp[1]$/, @{$phylums_u{$tmp[0]}}){
        push @{$phylums_u{$tmp[0]}}, $tmp[1];
        $uniq{$phylums{$tmp[1]}}++;
    }
    unless (grep /^$tmp[2]$/, @{$phylums_u{$tmp[0]}}){
        push @{$phylums_u{$tmp[0]}}, $tmp[2];
		$uniq{$phylums{$tmp[2]}}++;
    }
    my $line = join("\t", $phylums{$tmp[1]},$uniq{$phylums{$tmp[1]}},$uniq{$phylums{$tmp[1]}},$phylums{$tmp[2]},$uniq{$phylums{$tmp[2]}},$uniq{$phylums{$tmp[2]}});
    print LINK_OUT "$line\n";
}
close LINK;
close LINK_OUT;