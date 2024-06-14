use strict;

my ($in, $out) = @ARGV;
#perl $0 <input gbk path file> <output gbk info file>

open IN, $in or die $!;
open OUT, ">$out";
print OUT "Sample	GBK	GBK_size(Kb)	GBK_type	GBK_subtype	Leader_seq	Core_seq\n";
while(<IN>){
	chomp;
	my ($sample,$gbk,$path) = split("\t");
	my $len = 0;
	my $flag = 0;
	my (@types,@subs);
	open GBK, $path or die $!;
	while(<GBK>){
		chomp;
		if ($. == 1){
			my $ori_len = (split(" "))[2];
			$len = $ori_len/1000;
		}
		$flag = 1 if / protocluster /;
		$flag = 0 if / proto_core /;
		if ($flag && $_=~/category=/){
			my $type = ($_=~ /category="(.*)"/)[0];
			push @types, $type;
		}
		if ($flag && $_ =~ /product=/){
			my $sub = ($_=~ /product="(.*)"/)[0];
			push @subs, $sub;
		}
	}
	close GBK;
	my $t = join("/", @types);
	my $s = join("/", @subs);

	my ($l,$c) = qw/- -/;
	if ($t =~ /RiPP/){
		my (@cores, @leaders);
		open GBK, $path or die $!;
		my @lines = <GBK>;
		close GBK;
		for (my $i = 0; $i < scalar(@lines); $i++){
			if ($lines[$i] =~ /core_sequence=/){
				if ($lines[$i] =~ /.*"$/){
					my $seq = ($lines[$i]=~/\/core_sequence="(.*)"/)[0];
					push @cores, $seq;					
				}else{
					my $line = ($lines[$i]=~/\/core_sequence="(.*)/)[0];
					my $next_line = $lines[$i+1];
					$next_line =~s/[\s"]+//g;
					my $seq = $line . $next_line;
					push @cores, $seq					
				}
			}
			if ($lines[$i] =~ /leader_sequence=/){
				if ($lines[$i] =~ /.*"$/){
					my $seq = ($lines[$i]=~/\/leader_sequence="(.*)"/)[0];
					push @leaders, $seq;					
				}else{
					my $line = ($lines[$i]=~/\/leader_sequence="(.*)/)[0];
					my $next_line = $lines[$i+1];
					$next_line =~s/[\s"]+//g;
					my $seq = $line . $next_line;
					push @leaders, $seq					
				}
			}
		}
		$c = join(",", @cores) if scalar(@cores)!=0;
		$l = join(",", @leaders) if scalar(@leaders)!=0;
	}
	my $line = join("\t", $sample,$gbk,$len,$t,$s,$l,$c);
	print OUT "$line\n";
}
close IN;
close OUT;