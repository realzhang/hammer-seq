#!/usr/bin/env perl
# This script is used to extract methylation events from deduplicated mapped SAM files of Hammer-seq
# zhangzhuqiang@ibp.ac.cn
# last update: 2020.07.05
use strict;
use v5.10;

while(<>){
    my $r1 = $_;
    my $r2 = <>;
    my $s1 = ( split /:/, ( split /\s+/, $r1 )[-3] )[-1]; # z/h/x string
    my $s2 = ( split /:/, ( split /\s+/, $r2 )[-3] )[-1];
    my $flag =  ( split /\s+/, $r1 )[1];
    #say $flag;exit;
    #my $len = length($s1);

    # remove base of CpG that is half covered, ie., only C on one strand is covered. 
    # Allways happens at the start/end of the read. Four possibilities.
    if( uc( substr($s1, 0, 1) ) eq 'Z'    && uc( substr($s2, 1, 1) ) ne 'Z'   ){
        substr($s1, 0, 1) = '.';
    }

    if( uc( substr($s1, -1, 1) ) eq 'Z'    && uc( substr($s2, -2, 1) ) ne 'Z'   ){
        substr($s1, -1, 1) = '.';
    }

    if( uc( substr($s2, 0, 1) ) eq 'Z'    && uc( substr($s1, 1, 1) ) ne 'Z'   ){
        substr($s2, 0, 1) = '.';
    }

    if( uc( substr($s2, -1, 1) ) eq 'Z'    && uc( substr($s1, -2, 1) ) ne 'Z'   ){
        substr($s2, -1, 1) = '.';
    }

    my ($chr, $pos) = (split /\s+/, $r1)[2,3];
    my $zZ1 = ( $s1 =~ tr/zZ/zZ/ );
    my $zZ2 = ( $s2 =~ tr/zZ/zZ/ );
    next if ( $zZ1==0 || $zZ1 != $zZ2); # skip reads with no CpG
    my $ss1 =  $s1; my $ss2 = $s2;
    $ss1 =~ s/h|x|u|\.//ig ; # only keep Z/z
    $ss2 =~ s/h|x|u|\.//ig ;
    my @pos1; my @pos2;
    while( $s1 =~ /z/ig ){
        push @pos1, pos($s1);
    }
    while( $s2 =~ /z/ig ){
        push @pos2, pos($s2);
    }
    #say join("\t", "pos2:", @pos2);

    my $Z1 = ( $ss1 =~ tr/Z/Z/ );
    my $Z2 = ( $ss2 =~ tr/Z/Z/ );

    my $pdx = 0; # id of pd1 / pd2, upper arm or bottom arm of the replication fork; See illustration in Ming et al. Cell Res 2020.
    my $Sp; my $Sd; # Z/z/X/x/H/h string for parents and daughter strand;

    if($Z1 > $Z2){ # R1 has higher methylation
        $Sp = $ss1; $Sd = $ss2;
        #say join("\t", $chr, $pos, $pos+$len, $ss1, $ss2, $Z1, $z1, $Z2, $z2, $zZ1);
        if($flag == 0){ # + strand
            # R1 -> Pw; R2->Dc
            #print $Pw $r1;
            #print $Dc $r2;
            $pdx = 1;
        } else { # - strand
            # R1 -> Pc; R2 -> Dw
            #print $Pc $r1;
            #print $Dw $r2;
            $pdx = 2;
        }
    }

    elsif($Z1 < $Z2){ # R2 has higher methylation
       $Sp = $ss2; $Sd = $ss1;
        #say join("\t", $chr, $pos, $pos+$len, $ss2, $ss1, $Z2, $z2, $Z1, $z1, $zZ1);
        if($flag == 0){ # + strand        
            # R1 -> Dw; R2 -> Pc
            #print $Dw $r1;
            #print $Pc $r2;
            $pdx = 2;
        } else { # - strand
            # R1 -> Dc; R2 -> Pw
            #print $Dc $r1;
            #print $Pw $r2;
            $pdx = 1;
        }
    }

    else{ # Z1 == Z2, R1/R2 have the same methylation
        if($flag == 0){ # + strand
            # R1 -> Wat; R2 -> Crick
            my $ii = int(rand(2)); # classify randomly
            if( $ii == 0 ){
                #print $Pwe $r1; print $Dce $r2;
                $pdx = 1;
                $Sp = $ss1; $Sd = $ss2;
            }
            else{
                #print $Dwe $r1; print $Pce $r2;
                $pdx = 2;
                $Sp = $ss2; $Sd = $ss1;
            }
        } else { # - strand
            # R1 -> Crick; R2 -> Wat;
            my $ii = int(rand(2));
            if( $ii == 0 ){
                #print $Pce $r1; print $Dwe $r2;
                $pdx = 2;
                $Sp = $ss1; $Sd = $ss2;
            }
            else{
                #print $Dce $r1; print $Pwe $r2;
                $pdx = 1;
                $Sp = $ss2; $Sd = $ss1;
            }
        }
    }

    for my $i ( 0 .. $#pos1 ){
        my $offset = ( $pos1[$i] < $pos2[$i] ? $pos1[$i] : $pos2[$i] ) - 1;
        say join( "\t", $chr, $pos+$offset, $pdx, (split '', $Sp)[$i], (split '',$Sd)[$i] );
    }

}
# end of the file.
