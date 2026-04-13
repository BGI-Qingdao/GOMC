use strict;
# open IN, "trees.txt" or die $!;
open IN, $ARGV[0] or die $!;
my (%nodes, %bgcs, @gcfs, %uniq);
my $count = -1;
while(<IN>){
    chomp;
    $_ =~ s/\s//g;
    $_ =~ s/\(//g;
    $_ =~ s/\)//g;
    my ($node,$dist) = split(":", $_);
    $nodes{$count} = $node;

    my ($left,$right) = split(",", $node);
    get_bgcs($left, $count);
    get_bgcs($right, $count);
    # my $line = join(",", @{$bgcs{$count}});
    # print "$count\t$line\n";

    if ($dist>0.2){
        $uniq{$count} = 1;
        unless (exists $uniq{$left}){
            push @gcfs, $left;
        }
        unless (exists $uniq{$right}){
            push @gcfs, $right;
        }
    }
    $count += -1;
}
close IN;
# my $line = join(",", @gcfs);
# print "$line\n";

# open OUT, ">BGC_membership.txt";
open OUT, ">$ARGV[1]";
print OUT "\tgcf\n";
my $gcf = 1;
my $bgc;
# foreach my $node_site (0..$#gcfs){
foreach my $node (@gcfs){
    if ($node>=0){
        $bgc = $node+1;
        print OUT "$bgc\t$gcf\n";
    }else{
        foreach my $sbgc (@{$bgcs{$node}}){
            $bgc = $sbgc+1;
            print OUT "$bgc\t$gcf\n";
        }
    }
    $gcf++;
}
close OUT;


sub get_bgcs{
    my ($in, $count) = @_;
    if ($in>=0){
        push @{$bgcs{$count}}, $in;
    }else{
        push @{$bgcs{$count}}, @{$bgcs{$in}}
    }
    return @{$bgcs{$count}};
}