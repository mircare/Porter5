#!/usr/bin/perl

# Transforms a BLASTP output into a flat MA file containing in the first line
# the number of alignments and the sequences in the following lines
#


$filen=$ARGV[0];
$fileo=$ARGV[1];
$fasta=$ARGV[2];

open(fast,"<$fasta");
@testo_temp=<fast>;
shift(@testo_temp);
while (@testo_temp) {
	chomp(@testo_temp);
	$truequery.=shift(@testo_temp);
}
chomp($truequery);
close(fast);


#open(fi,"<$filen");
#@testo=<fi>;
#close(fi);

$conta=0;
$contaq=0;
$contas=0;

$query="";
$sbjct="";
$aperto=0;
$first=1;

$offset=0;
$l=0;

$ofile="";

#while (@testo) {
#	$linea=shift(@testo);

$temp = $fileo.'.app';
open(o,">$temp");
print o "$truequery\n";
$l=length($truequery);
#print "$l\n";
close(o);

open(fi,"<$filen");
@text = <fi>;
close fi;


if ($#text<0) {
	open(o,">$fileo");
	print o "1\n$truequery";
	close(o);
	exit;
}


for ($pp=0;$pp<=$#text;$pp++) {

	$linea = $text[$pp];
	chomp($linea);

#	if ($conta-$contaq==5 && $aperto && $first) {
#		print o "$query\n";
#		$l=length($query);
#		$first=0;
#		}
	if ($conta-$contas==5 && $aperto) {

		while (index($query,"-")!=-1) {
			$pos=index($query,"-");
			$query=substr($query,0,$pos).substr($query,$pos+1);
			$sbjct=substr($sbjct,0,$pos).substr($sbjct,$pos+1);
			}
		$temp=(".") x $offset;
		$sbjct = $temp.$sbjct;
		if (length($sbjct)>$l) {$sbjct = substr($sbjct,0,$l);}
		$temp=(".") x ($l-length($sbjct));
		$sbjct = $sbjct.$temp;
		while (index($sbjct,"-")!=-1) {
			$pos=index($sbjct,"-");
			$sbjct=substr($sbjct,0,$pos).'.'.substr($sbjct,$pos+1);
			}

		print o "$sbjct\n";
		$query="";
		$sbjct="";
		}


	if ($linea =~ /^Query= ([\000-\256]*)/) {
		$temp=$1;
		if ($temp =~ /^([\000-\256]*)([\.]+)([\000-\256]*)/) {
			$temp=$1;
			}
		$ofile=$fileo;
#		print "$ofile ";
		$temp = $fileo.'.app';
		open(o,">>$temp");
		$aperto=1;
		$first=1;
		}
	elsif ($pp == $#text) { #($linea =~ /^S2: ([\000-\256]*)/) {
		close(o);
		$tmp = $fileo.'.app';
		open(fa,"<$tmp");
		@testo=<fa>;
		close(fa);
		open(o,">$ofile");
		$num=$#testo+1;
		print o "$num\n";
		while (@testo) {
			$app=shift(@testo);
			print o $app;
			}
		close(o);
		$aperto=0;
		}
	elsif ($linea =~ /^Query  ([0-9]+) ([\000-\256]*)/) {
		if ($conta-$contaq!=4) 
			{$offset=$1-1;}
#		print $offset," ";
		$temp=$2;
		if ($temp =~ /^([ ]*)([A-Z\-]+)  ([0-9]+)/) {
			$temp=$2;
			}
#		print o $temp,"\n";
		if ($conta-$contaq==4) {
			$query .= $temp;
			}
		else {$query=$temp;}
		$contaq=$conta;
		}
	elsif ($linea =~ /^Sbjct  ([\000-\256]*)/) {
		$temp=$1;
		if ($temp =~ /^([\000-\256]*)([ ]+)([\000-\256]*)  ([\000-\256]*)/) {
			$temp=$3;
			}
#		print o $temp,"\n";
		if ($conta-$contas==4) {
			$sbjct .= $temp;
			}
		else {$sbjct=$temp;}
		$contas=$conta;
		}
	$conta++;
	}

#close(o);

system("rm $tmp");
close(fi);




