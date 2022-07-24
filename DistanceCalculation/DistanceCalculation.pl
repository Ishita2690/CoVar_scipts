#! /usr/bin/perl

use List::Util qw( min max );

$struc1=@ARGV[0];
$pdb=@ARGV[1];
$cha=@ARGV[2];
$cha2=@ARGV[3];
$Thresh=@ARGV[4];

@bon=split(/=/,$cha);
@bon2=split(/=/,$cha2);
chomp $bon[0];
chomp $bon2[0];
chomp $bon[1];
chomp $bon2[1];
$dat{$bon2[0]}="$bon2[1]";
$dat{$bon[0]}="$bon[1]";

#@st1=`cat $struc1`;
open $str11, "$struc1" or die "Error: reading position file";
@st1=(<$str11>);
#@pdg = `cat $pdb`;

open $PDB, "$pdb" or die "Error: reading structure file";
@pdg =(<$PDB>);

@mon=split(/\t/,$st1[0]);
chomp $mon[0];
@von=split(/\,/,$mon[0]);
chomp $von[0];
chomp $von[1];

@dem=split(/\[/,$von[0]);
chomp $dem[1];
@nem=split(/\]/,$dem[1]);
chomp $nem[0];
@dem2=split(/\[/,$von[1]);
chomp $dem2[1];
@nem2=split(/\]/,$dem2[1]);
chomp $nem2[0];

$nem[0]=~s/\s+//g;
$nem2[0]=~s/\s+//g;

for($l=1;$l<=$#st1;$l++){
	 chomp $st1[$l];
	 @bos=split(/\t/,$st1[$l]);
	 chomp $bos[0];
	 chomp $bos[1];
	 chomp $bos[2];
	 chomp $bos[3];
	 chomp $bos[5];
	 @pos=split(/\,/,$bos[0]);
	 chomp $pos[1];
	 chomp $pos[0];
	 $k=0;
	while($k<=$#pdg)
	{
		 chomp $pdg[$k];
		 @wer=split(/\s+/,$pdg[$k]);
		 chomp $wer[4];
		 chomp $wer[0];
		 chomp $wer[5];

##checkpoint to match PDB chain and residue information
		if( $nem[0] eq $bon[0] )
		{
			 if( $wer[0] eq "ATOM" && $wer[4] eq $dat{$bon[0]} && $wer[5] eq $pos[0] )
			{
				push( @pdba, "$pdg[$k]");
			}
			if( $wer[0] eq "ATOM" && $wer[4] eq $dat{$bon2[0]} && $wer[5] eq $pos[1] )
			{
	                        push (@pdbb, "$pdg[$k]");
			}
		}
		elsif( $nem2[0] eq $bon[0] )
		{
		if( $wer[0] eq "ATOM" && $wer[4] eq $dat{$bon2[0]} && $wer[5] eq $pos[0] )
			{
                        push( @pdba, "$pdg[$k]");
                	}
                if( $wer[0] eq "ATOM" && $wer[4] eq $dat{$bon[0]} && $wer[5] eq $pos[1] )
			{
                        push (@pdbb, "$pdg[$k]");

	
			}
		}

	$k++;
	}

		


    $m=0;
    while($m<=$#pdba)
    {
     chomp $pdba[$m];
     @xx=split(/\s+/,$pdba[$m]);
     $x1=@xx[6];$y1=@xx[7];$z1=@xx[8];
     $n=0;
     	while($n<=$#pdbb)
     	{
      	chomp $pdbb[$n];
      	@yy=split(/\s+/,$pdbb[$n]);
      	$x2=@yy[6];$y2=@yy[7];$z2=@yy[8];
      	$x=($x1-$x2)**2;
      	$y=($y1-$y2)**2;
      	$z=($z1-$z2)**2;
      	$distance=sprintf("%0.2f",((sqrt($x+$y+$z))));
       	push(@diste,"$distance");
       	undef @yy;
        $n++;
     }
     undef @xx;
     $m++;
  }


	
	$mind= min @diste;
	$MinD{$pos[0]}{$pos[1]}="$mind";
	$dog{$pos[0]}{$pos[1]}="$bos[1]";
	$cos{$pos[0]}{$pos[1]}="$bos[2]";
	$zs{$pos[0]}{$pos[1]}="$bos[3]";
	$Sres{$pos[0]}{$pos[1]}="$bos[5]";
	keys %MinD;
	keys %cos;
	keys %zs;
	keys %dog;
	keys %Sres;
	undef @pdba;
	undef @pdbb;
	undef @diste;
}


open(kol,">$struc1-Min_Dist");

print kol"Position 1 [$nem[0]],Position 2 [$nem2[0]] (Reference structure)\tPosition 1 [$nem[0]],Position 2 [$nem2[0]] (Reference sequence)\tMinimum Distance Between Pairs\tCo-Var Score\tZ-Score\tResidue 1, Residue 2\n";


open(DoInt,">$struc1-Min_Dist_Thresh");

print DoInt"Position 1 [$nem[0]],Position 2 [$nem2[0]] (Reference structure)\tPosition 1 [$nem[0]],Position 2 [$nem2[0]] (Reference sequence)\tMinimum Distance Between Pairs\tCo-Var Score\tZ-Score\tResidue 1, Residue 2\n";



foreach $we (keys %MinD)
{

	foreach $re (keys %{$MinD{$we}})
	{
		chomp $we;
		chomp $re;

		print kol"$we,$re\t$dog{$we}{$re}\t$MinD{$we}{$re}\t$cos{$we}{$re}\t$zs{$we}{$re}\t$Sres{$we}{$re}\n";

		if($MinD{$we}{$re} <= $Thresh && $MinD{$we}{$re} != "")
		{
		print DoInt"$we,$re\t$dog{$we}{$re}\t$MinD{$we}{$re}\t$cos{$we}{$re}\t$zs{$we}{$re}\t$Sres{$we}{$re}\n";

		}
		else
		{
		next;
		}


	}

}
undef @st1;
undef @comb;
undef @pdg;
undef @pos;
undef @wer;
undef %MinD;
undef %dog;
undef %zs;
undef %cos;
undef %Sres;
close $PDB;
close $Str1;
close(kol);
close(DoInt);
exit;
