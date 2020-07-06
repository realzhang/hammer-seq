#!/usr/bin/env perl
use v5.10;
use strict;

my %stat;
my @stat_overall;

while(<>){
    chomp;
    my @dd = split /\s+/, $_;
    my ($chr,$start, $end, $ss) = @dd[0,1,2,  $#dd - 1];
    my ($ZZ1, $Zz1, $zz1, $zZ1, $ZZ2, $Zz2, $zz2, $zZ2) = split /,/, $ss;

    my $k = join(".", $chr, $start, $end);

    $stat_overall[0] = ($stat_overall[0]//0) + $ZZ1;
    $stat_overall[1] = ($stat_overall[1]//0) + $Zz1;
    $stat_overall[2] = ($stat_overall[2]//0) + $zz1;
    $stat_overall[3] = ($stat_overall[3]//0) + $zZ1;
    $stat_overall[4] = ($stat_overall[4]//0) + $ZZ2;
    $stat_overall[5] = ($stat_overall[5]//0) + $Zz2;
    $stat_overall[6] = ($stat_overall[6]//0) + $zz2;
    $stat_overall[7] = ($stat_overall[7]//0) + $zZ2;

    if( exists $stat{$k} ){
        $stat{$k}[0] += $ZZ1;
        $stat{$k}[1] += $Zz1;
        $stat{$k}[2] += $zz1;
        $stat{$k}[3] += $zZ1;

        $stat{$k}[4] += $ZZ2;
        $stat{$k}[5] += $Zz2;
        $stat{$k}[6] += $zz2;
        $stat{$k}[7] += $zZ2;
    }
    else {
        $stat{$k} = [$ZZ1, $Zz1, $zz1, $zZ1, $ZZ2, $Zz2, $zz2, $zZ2];
    }
}

# Print the header of ouput in file 'events_and_maintenance_by_ROI.txt'
say "Chr.start.end\tUp.MM\tUp.MU\tUp.UU\tUp.UM\tLo.MM\tLo.MU\tLo.UU\tLo.UM\tMaint.ratio\tDenovo.ratio";
for ( keys %stat ) {
    my @counts = @{$stat{$_}};
    my $mat = ($counts[0]+$counts[1]+$counts[4]+$counts[5])>=5 ? ($counts[0] + $counts[4]) / ($counts[0]+$counts[1]+$counts[4]+$counts[5]) : 'na';
    my $den = ($counts[2]+$counts[3]+$counts[6]+$counts[7])>=5 ? ($counts[3] + $counts[7]) / ($counts[2]+$counts[3]+$counts[6]+$counts[7]) : 'na';
    say join("\t", $_, @counts, $mat, $den );
}

my $mat_overall = ($stat_overall[0] + $stat_overall[4]) / ($stat_overall[0]+$stat_overall[1]+$stat_overall[4]+$stat_overall[5]);
my $den_overall = ($stat_overall[3] + $stat_overall[7]) / ($stat_overall[2]+$stat_overall[3]+$stat_overall[6]+$stat_overall[7]);

say STDERR "For your ROIs:\nThe overall maintenance ratio is:\t$mat_overall\nThe overall de novo methylation ratio is:\t$den_overall.";
