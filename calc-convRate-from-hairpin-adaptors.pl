#!/usr/bin/env perl
# This script is used for Hammer-seq analysis, calculate BS conversion rate from unmethylated hairpin adaptors.
# zhangzhuqiang@ibp.ac.cn
# last update: 2020.07.03;
use strict;
use v5.10;

my $reads_num = 0;
my $numC = 0;

while(<>){
    # get seq from the 2nd of every 4 lines;
    my $seq = <>;
    <>;<>; # skip 3rd, 4th lines.
    chomp $seq;

    if ($seq =~ /([CT]GGATGAG[CT]AAGTGAAG[CT]T[CT]AT[CT][CT]GAT[CT])/){
        my $m = $1;
        $reads_num += 1;
        my $C = $m =~ tr/C/C/;
        $numC += $C;
    }
    print  "$. lines.\n" if ($. % 10000000 == 0);
}

# output to STDOUT
print "\n";
say "Total hairpins found:\t $reads_num";
say "Total C in hairpins:\t", $reads_num*7;
say "Non converted C:\t $numC";
say "Non converted C (%):\t", sprintf("%.2f", $numC / $reads_num / 7 * 100), '%';
say "Convertion rate (%):\t", sprintf("%.2f", (1-$numC / $reads_num / 7)*100), '%';

# output to file named 'conv_rate'
open my $FH, ">>conv_rate";
say $FH "Total hairpins found:\t $reads_num";
say $FH "Total C in hairpins:\t", $reads_num*7;
say $FH "Non converted C:\t $numC";
say $FH "Non converted C (%):\t", sprintf("%.2f", $numC / $reads_num / 7 * 100), '%';
say $FH "Convertion rate (%):\t", sprintf(("%.2f", 1-$numC / $reads_num / 7)*100), '%';

# End of the file.
