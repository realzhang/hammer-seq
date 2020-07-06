#!/usr/bin/env perl
# This script is used to count frequencies of methylation events on each genomic coordinates.
use strict;
use v5.10;

my %d;
while(<>){
    chomp;
    my ($chr, $pos, $pdx, $P, $D) = split /\t/,$_;
    #say "chr is $chr, pos is $pos";exit;
    $d{"$chr.$pos"}{$pdx}{"$P$D"} = $d{"$chr.$pos"}{$pdx}{"$P$D"} ? $d{"$chr.$pos"}{$pdx}{"$P$D"}+1 : 1;
}

say STDERR "Now print:";
open my $O_bed, ">>events_count.bed" or die $!;
open my $O_txt, ">>events_count.txt" or die $!;

# print header to txt file
say $O_txt "chr\tstart\tend\tUp.MM\tUp.MU\tUp.UU\tUp.UM\tLo.MM\tLo.MU\tLo.UU\tLo.UM";

for my $k ( keys %d ){
    my ($chr, $pos) = split(/\./,$k);

    say $O_bed join("\t", $chr, $pos, $pos+1,
        join(',',    
        $d{$k}{1}{'ZZ'} // 0,
        $d{$k}{1}{'Zz'} // 0,
        $d{$k}{1}{'zz'} // 0,
        $d{$k}{1}{'zZ'} // 0,
        $d{$k}{2}{'ZZ'} // 0,
        $d{$k}{2}{'Zz'} // 0,
        $d{$k}{2}{'zz'} // 0,
        $d{$k}{2}{'zZ'} // 0 )
    );

    say $O_txt join("\t", $chr, $pos, $pos+1,
        $d{$k}{1}{'ZZ'} // 0,
        $d{$k}{1}{'Zz'} // 0,
        $d{$k}{1}{'zz'} // 0,
        $d{$k}{1}{'zZ'} // 0,
        $d{$k}{2}{'ZZ'} // 0,
        $d{$k}{2}{'Zz'} // 0,
        $d{$k}{2}{'zz'} // 0,
        $d{$k}{2}{'zZ'} // 0 );
}
