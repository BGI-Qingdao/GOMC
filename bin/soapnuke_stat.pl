#!/usr/bin/perl

my ($basic,$filter)=@ARGV;

my %hash=();
open IN,$basic or die $!;
while (<IN>)
{
	chomp;
	my @arr=split /\t/,$_;
	next if($arr[0] =~ /^Item/);
	my %num=();
	my %fre=();
	foreach my $i (1..4)
	{
		if($arr[$i] =~ /^(\d+)\s+\(\s*([\d\.]+)\%\)$/)
		{
			$num{$i}=$1;
			$fre{$i}=$2;
		}
	}
	if($arr[0] =~ /^Total number of reads/)
	{
		$hash{'tr_r'}=$num{1}+$num{3};
		$hash{'tr_c'}=$num{2}+$num{4};
	}
	if($arr[0] =~ /^Total number of bases/)
	{
		$hash{'tb_r'}=$num{1}+$num{3};
		$hash{'tb_c'}=$num{2}+$num{4};
	}
	if($arr[0] =~ /^Number of base N/)
	{
		$hash{'n1_r'}=$num{1};
		$hash{'n1_c'}=$num{2};
		$hash{'n2_r'}=$num{3};
		$hash{'n2_c'}=$num{4};
	}
	if($arr[0] =~ /^Number of base C/)
	{
		$hash{'c1_r'}=$fre{1};
		$hash{'c1_c'}=$fre{2};
		$hash{'c2_r'}=$fre{3};
		$hash{'c2_c'}=$fre{4};
	}
	if($arr[0] =~ /^Number of base G/)
	{
		$hash{'g1_r'}=$fre{1};
		$hash{'g1_c'}=$fre{2};
		$hash{'g2_r'}=$fre{3};
		$hash{'g2_c'}=$fre{4};
	}
	if($arr[0] =~ /^Number of base calls with quality value of 20 or higher/)
	{
		$hash{'q20_1r'}=$fre{1};
		$hash{'q20_1c'}=$fre{2};
		$hash{'q20_2r'}=$fre{3};
		$hash{'q20_2c'}=$fre{4};
	}
	if($arr[0] =~ /^Number of base calls with quality value of 30 or higher/)
	{
		$hash{'q30_1r'}=$fre{1};
		$hash{'q30_1c'}=$fre{2};
		$hash{'q30_2r'}=$fre{3};
		$hash{'q30_2c'}=$fre{4};
	}
}
close IN;

open IN,$filter or die $!;
while (<IN>)
{
	chomp;
	my @arr=split /\t/,$_;
	next if($arr[0] eq 'Item');
	if($arr[0] =~ /^Read with n rate exceed:/)
	{
		$hash{'nr'}=$arr[1];
	}
	if($arr[0] =~ /^Reads with low quality/)
	{
		$hash{'lr'}=$arr[1];
	}
	if($arr[0] =~ /^Reads with adapter/)
	{
		$hash{'ar'}=$arr[1];
	}
}
close IN;

my $gc_1r=$hash{'g1_r'}+$hash{'c1_r'};
my $gc_1c=$hash{'g1_c'}+$hash{'c1_c'};
my $gc_2r=$hash{'g2_r'}+$hash{'c2_r'};
my $gc_2c=$hash{'g2_c'}+$hash{'c2_c'};

print "Type\tRaw data\tClean data\n";
print "Number of Reads\t$hash{'tr_r'}\t$hash{'tr_c'}\n";
print "Data Size\t$hash{'tb_r'}\t$hash{'tb_c'}\n";
print "N of fq1\t$hash{'n1_r'}\t$hash{'n1_c'}\n";
print "N of fq2\t$hash{'n2_r'}\t$hash{'n2_c'}\n";
print "GC(%) of fq1\t$gc_1r\t$gc_1c\n";
print "GC(%) of fq2\t$gc_2r\t$gc_2c\n";
print "Q20(%) of fq1\t$hash{'q20_1r'}\t$hash{'q20_1c'}\n";
print "Q20(%) of fq2\t$hash{'q20_2r'}\t$hash{'q20_2c'}\n";
print "Q30(%) of fq1\t$hash{'q30_1r'}\t$hash{'q30_1c'}\n";
print "Q30(%) of fq2\t$hash{'q30_2r'}\t$hash{'q30_2c'}\n";
print "Discard Reads related to N\t$hash{'nr'}\n";
print "Discard Reads related to low qual\t$hash{'lr'}\n";
print "Discard Reads related to Adapter\t$hash{'ar'}\n";
