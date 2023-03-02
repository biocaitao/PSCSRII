use strict;

my $binlen=10;
if(@ARGV<2){
  die "perl filter_miRNA.pl gencode.v32.annotation.gtf  seq.map;";
}

my $of = $ARGV[1] . ".deg";
open(OF, ">$of") or die "$!";
my %infor = ();
my $pos = "";
my $pid = "";
open(INFL, "samtools view $ARGV[1] |" ) or die "$!";
while(<INFL>){
  my @t =split /\t/;
  my @t2 = split /_/, $t[0];
  my $chr = $t[2];
  my $strand = "+";
  if($t[1] == 16 or $t[1] == 272){
    $strand = "-";
  }
  
  my $cell = $t2[1];
  my $umi = $t2[2];
# my @t3 =split / /, $t[-1];
#consider pos
#for each moleculor, for each cell, for each pos,  for each UMI
  if($pos == $t[3]){
    $pid = $t[3] . $t[5];
    if(! exists $infor{$pid}{$strand}{$cell}{$umi}){
     print OF $_;
     $infor{$pid}{$strand}{$cell}{$umi} = 1;
     my @misumis = mismatchstr($umi);
     $infor{$pid}{$strand}{$cell}{$_}++ for (@misumis);
    }
  }else{
    $pos = $t[3];
    $pid = $t[3] . $t[5];
    %infor=();
    $infor{$pid}{$strand}{$cell}{$umi} = 1;
    my @misumis = mismatchstr($umi);
    $infor{$pid}{$strand}{$cell}{$_}++ for (@misumis);
    print OF $_;
  }

}
close(INFL);
close(OF);


my %ncRNA=();
open(INFL, $ARGV[0]) or die "$!";
while(<INFL>){
  chomp;
  my @t =split /\t/;
  my ($id) = $t[2] ; 
#-1 shift to 0 base
  my $bin_start = $t[3]-$binlen-1;
  my $bin_end = $t[4] -1;
  for(my $i=$bin_start; $i<=$bin_end; $i++){
    $ncRNA{$t[0]}{$t[6]}{$i}{$id} = 1;
  }

}
close(INFL);


open(INFL,  $of) or die "$!";
while(<INFL>){
  chomp;
  my @t=split /\t/;
#  next if($t[11] ne "NH:i:1");
  my $po = $t[3];
  my $strand = "+";
  if($t[1] == 16 or $t[1] == 272){
    $strand = "-";
  }
  if(exists $ncRNA{$t[2]}{$strand}{$po}){
       print join("\t", @t);
       print "\t", $t[11], "\t";
       print join(" ", keys %{$ncRNA{$t[2]}{$strand}{$po}}), "\n";

  }
 
}
close(INFL);


sub mismatchstr {
my ($test) =@_;
my @str =split //, $test;
my @str2 =();
my @dict = ("A", "G", "C", "T", "N");;
for(my $i=0; $i<@str; $i++){
  my @t = grep {$str[$i] ne $_} @dict;
  foreach my $e (@t){

    my @s = @str;
    $s[$i] = $e;
    my $s2 = join("", @s);
    push @str2, $s2;
  }

}
return @str2;
}

