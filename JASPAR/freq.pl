#!/usr/bin/perl

while(my $file =<MA0*>)
{
chomp($file);
open(IN,$file);
open(OUT,">>all.mat");
print OUT ">".substr($file,0,index($file,"."))."\n";
my @nuc=();
$nuc[0]=();
$nuc[1]=();
$nuc[2]=();
$nuc[3]=();
my $c=0;
while (my $line=<IN>)
{
	chomp($line);
	$nuc[$c++]=($line=~/(\d+)/g);
}
	my $sum = $aga[0]+$aga[1]+$aga[2]+$aga[3];
	printf OUT ("%f %f %f %f\n",$aga[0]/$sum,$aga[1]/$sum,$aga[2]/$sum,$aga[3]/$sum);
print OUT "<\n";
close(OUT);
close(IN);
}
