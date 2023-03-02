use strict;

if(@ARGV<3){
  die "perl stat_reads.pl processed.fa processed2.fa sRNAAlignment.bam";

}
#reads with adapters
my %count1 =();
open(INFL, $ARGV[0]) or die "$!";
while(<INFL>){
  if(/^\@/){
    chomp;
    my @t =split /_/, $_;
    $count1{$t[1]}{$t[0]} =1;
  }
}
close(INFL);

#reads with correct length
my %count2 =();
open(INFL, $ARGV[1]) or die "$!";
while(<INFL>){
  if(/^\@/){
    chomp;
    my @t =split /_/, $_;
    $count2{$t[1]}{$t[0]} =1;
  }
}
close(INFL);

#matched small reads counts
my %count3 =();
#matched UMI counts
my %count4 =();
open(INFL,  "samtools view $ARGV[2] |") or die "$!";
while(<INFL>){
  my @t=split /\t/;
  if(($t[1] == 0) or ($t[1] == 16)){
    my @t2 =split /_/, $t[0];
    $count3{$t2[1]}{$t2[0]} =1;
    $count4{$t2[1]}{$t2[2]} =1;
  }
}
close(INFL);


foreach my $cell (keys %count1){
  print $cell, "\t";
  print scalar keys %{$count1{$cell}} || 0;
  print "\t";
  print scalar keys %{$count2{$cell}} || 0; 
  print "\t";
  print scalar keys %{$count3{$cell}} || 0;
  print "\t";
  print scalar keys %{$count4{$cell}} || 0;
  print "\n";

}
