#!/dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/envs/bigslice/bin/perl
use strict;
use Statistics::Descriptive;

open OUT, ">Variance_Info.txt";
print OUT "Level\tName\tGCF numbers of lower-ranked taxa\tvariance\n";
# open OUT, ">variance.result";
# print OUT "Type\tvariance\n";
# foreach my $file (glob "/dellfsqd2/ST_OCEAN/USER/zhouchanghao/metagenome/get_gcf/20220905/Output/Variance/result/*.txt"){
foreach my $file (glob "./result/*.txt"){
    # print "$file\n";
    my $tag = ($file =~ /.*\/[0-9]\.([a-z])_.*/)[0];
    $tag = trans($tag);
    my %info;
    open FILE, $file or die $!;
    while(<FILE>){
        chomp;
        my @tmp = split("\t");
        next if $tmp[2] =~ /_$/;
        $info{$tmp[0]}{$tmp[2]} = $tmp[3];
    }
    close FILE;

    my %need;
    foreach my $u (sort keys %info){
        if ((keys %{$info{$u}})>1){
            foreach my $d (sort keys %{$info{$u}}){
                push @{$need{$u}}, $info{$u}{$d};
            }
        }
    }

    foreach my $u (sort keys %need){
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(@{$need{$u}});
        my $v = $stat->variance();
        # my $line = join("\t", $tag,$v);
        # print OUT "$line\n";
        my $set = join(",", @{$need{$u}});
        my $line = join("\t", $tag,$u,$set,$v);
        print OUT "$line\n";
    }
}

sub trans{
    my ($in) = @_;
    my $out;
    $out = "Phyla" if $in eq "p";
    $out = "Classes" if $in eq "c";
    $out = "Orders" if $in eq "o";
    $out = "Families" if $in eq "f";
    $out = "Genera" if $in eq "g";
    return $out;
}