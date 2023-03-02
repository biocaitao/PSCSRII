use strict;

my %wl =();
if(@ARGV<2) {
  die "perl parseseq.pl whitelist.txt seq.gz";
}

open(INFL, "$ARGV[0]") or die"need whitelist.txt file $!";
while(<INFL>){
  chomp;
  my @t =split /\t/;
#   print $t[0], "\n";
  $wl{$t[0]}=$t[0];
  my @t2 =split /,/, $t[1];
  foreach my $id (@t2){
#     print $id, "\n";
     $wl{$id}=$t[0];
  }
}
close(INFL);

open(INFL1, "zcat $ARGV[1] |") or die "$!";

open(INFL2, "zcat $ARGV[2] |") or die "$!";


my $find1 = "TGTCACG(T+)(G*)";
my $find2 = "AAAAAAAA";

my $li=0;
while(!eof(INFL1)){
  my $id1 =  <INFL1>;
  my $id2 =  <INFL2>;

  my $seq = <INFL1>;
  chomp($seq);
  my $seqe =<INFL2>;
  my $t1 = <INFL1>;
  my $t2 = <INFL2>;

  my $q1 = <INFL1>;
  chomp($q1);
  my $q2 = <INFL2>;
  $li++;
  $id1=~ s/ (.*)\n//;;
#  print $id, "\t", $seq;
  
  if($seq=~/$find1?/){
#     print $seq;
       my $p1 = $-[0];
       my $s1 = $+[0];

       my $oseq = "";
       my $oq = "";
       my $barcode1 = substr($seq, 0, 6);
       my $barcode2 = substr($seqe, 0, 6);
       my $revb2 = reverse  $barcode2;
       $revb2 =~ tr/ATGCatgc/TACGtacg/;
       my $bc = $barcode1 . $revb2;
       my $umi = substr($seqe, 30, 12);

       if( ($seq=~/$find2/ )) {
         my $s2 = $-[0]-2-$s1;
         if($s2>0){
         $oseq =  substr($seq, $s1, $s2); 
         $oq = substr($q1, $s1, $s2);
         }
       }else{

         $oseq = substr($seq, $s1);
         $oq = substr($q1, $s1);
       }
       

       if(exists $wl{$bc}){

           print $id1, "_";
           print $wl{$bc} ,"_",$umi,  "\n";
           print $oseq,  "\n";
           print "+\n";
           print $oq, "\n";

      }
      
 }


}

close(INFL1);
close(INFL2);
print STDERR "processed $li raw reads\n";

sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
