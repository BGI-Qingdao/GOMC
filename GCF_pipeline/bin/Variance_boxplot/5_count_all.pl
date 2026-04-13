use strict;

my %info;
my %count;
open IN, "./1_GCF_info_all.txt" or die $!;
while(<IN>){
    chomp;
    my @tmp  =split("\t");
    my @levels = (split(";",$tmp[2]));
    my ($d,$p,$c,$o,$f,$g,$s) = (split(";",$tmp[2]));
    unless (grep /^$tmp[0]$/, @{$info{$p}{$c}{$o}{$f}{$g}{$s}}){
        push @{$info{$p}{$c}{$o}{$f}{$g}{$s}}, $tmp[0];
    }
    my $cl = join(";", $p,$c);
    my $ol = join(";", $cl,$o);
    my $fl = join(";", $ol,$f);
    my $gl = join(";", $fl,$g);
    my $sl = join(";", $gl,$s);
    unless (grep /^$tmp[0]$/, @{$count{$p}}){
        push @{$count{$p}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$count{$cl}}){
        push @{$count{$cl}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$count{$ol}}){
        push @{$count{$ol}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$count{$fl}}){
        push @{$count{$fl}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$count{$gl}}){
        push @{$count{$gl}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$count{$sl}}){
        push @{$count{$sl}}, $tmp[0];
    }
}
close IN;

open OUT, ">./All_levels_Info.txt";
print OUT "Phyla\tClasses\tOrders\tFamilies\tGenus\tSpecies\n";
foreach my $p (sort keys %info){
    foreach my $c (sort keys $info{$p}){
        foreach my $o (sort keys $info{$p}{$c}){
            foreach my $f (sort keys $info{$p}{$c}{$o}){
                foreach my $g (sort keys $info{$p}{$c}{$o}{$f}){
                    foreach my $s (sort keys $info{$p}{$c}{$o}{$f}{$g}){
                        my $s_count = $#{$info{$p}{$c}{$o}{$f}{$g}{$s}}+1;
                        my $cl = join(";", $p,$c);
                        my $ol = join(";", $cl,$o);
                        my $fl = join(";", $ol,$f);
                        my $gl = join(";", $fl,$g);
                        my $sl = join(";", $gl,$s);
                        my $pnum = $#{$count{$p}}+1;
                        my $cnum = $#{$count{$cl}}+1;
                        my $onum = $#{$count{$ol}}+1;
                        my $fnum = $#{$count{$fl}}+1;
                        my $gnum = $#{$count{$gl}}+1;
                        my $snum = $#{$count{$sl}}+1;
                        my $sline = join("\t", "$p($pnum)","$c($cnum)","$o($onum)","$f($fnum)","$g($gnum)","$s($snum)");
                        print OUT "$sline\n";
                    }
                }
            }
        }
    }
}
close OUT;
