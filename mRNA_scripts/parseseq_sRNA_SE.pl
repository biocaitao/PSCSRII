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

open(INFL, "zcat $ARGV[1]|") or die "$!";
my $find1 = "TTG([ACGT][ACGT])TTG([ACGT][ACGT])TTG([ACGT][ACGT])";

my $find2 = "CTGTAG";

my $li=0;
while(<INFL>){
  my $id =  $_;
  my $seq = <INFL>;
  my $t = <INFL>;
  my $q = <INFL>;
  $li++;
  $id=~ s/ (.*)\n//;;
#  print $id, "\t", $seq;
  my $seq1 = substr($seq, 6, 18);
  my $seq2 = substr($seq, 7, 18);
  my $seq3 = substr($seq, 8, 18);
  my $seq4 = substr($seq, 9, 18);
  my $ref = "GCAGTGGTAGTTATCGCG";
  my $seq5 = "";
 
  my $s1 =0;
  
  if($seq=~/$find1/  and $-[0] >4 ){
#     print $seq, $s1;
       my $um1 = $1;
       my $um2 = $2;
       my $um3 = $3;
       $s1 =  $+[0];
       if( ($seq=~/$find2/ )and ($-[0] < length($seq) ) ){
         my $um4 = substr($seq, $-[0]-2, 2);
         my $s2 = $-[0]-2 - $s1;
#         print $-[0]-2, " ", $s1, "\n";
         if($s2 >=0){
         my $barcode1 = substr($seq, 0, 6);
         my $barcode2 = substr($seq, $-[0]+24, 6);
         my $bc = $barcode1 . $barcode2;

         if(exists $wl{$bc}){

           print ">", $id, "_";
           print $wl{$bc}, "_", $um1, $um2, $um3, $um4 ,"\n";
           print substr($seq, $s1, $s2), "\n";
         }
         }
      
     }

 }

}

close(INFL);
print STDERR "processed $li raw reads\n";

sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
