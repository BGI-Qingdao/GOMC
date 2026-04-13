use strict;

my $need = $ARGV[0];
my %info;
my %counts;
open IN, "./1_GCF_info_all.txt" or die $!;
while(<IN>){
    chomp;
    my @tmp  =split("\t");
    my $u = (split(";",$tmp[2]))[$need];
    my $d = (split(";",$tmp[2]))[$need+1];
    unless (grep /^$tmp[0]$/, @{$info{$u}{$d}}){
        push @{$info{$u}{$d}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$counts{$u}}){
        push @{$counts{$u}}, $tmp[0];
    }
}
close IN;


my $name = &outname($need);
`mkdir result` unless -e "result";
open OUT, ">result/$name";
foreach my $u (sort keys %info){
    my $u_num = $#{$counts{$u}}+1;
    foreach my $d (sort keys $info{$u}){
        my $num = $#{$info{$u}{$d}}+1;
        print OUT "$u\t$u_num\t$d\t$num\n";
    }
}
close OUT;

sub outname{
    my ($need) = @_;
    my $out;
    if ($need == 1){
        $out = "1.p_c_result.txt";
    }elsif($need == 2){
        $out = "2.c_o_result.txt";
    }elsif($need == 3){
        $out = "3.o_f_result.txt";
    }elsif($need == 4){
        $out = "4.f_g_result.txt";
    }elsif($need == 5){
        $out = "5.g_s_result.txt";
    }
    return $out;
}
