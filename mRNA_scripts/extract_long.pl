use strict;
open(INFL, $ARGV[0]) or die "perl extract_sRNA.pl processed.fa";
my $i=0;
my $j=0;
while(<INFL>){
  my $id = $_;
  my $seq = <INFL>;
  my $t = <INFL>;
  my $q = <INFL>;
  my $l = length($seq);
  $i++;
#  print $id, $seq, "\n";
  if( $l>=30 ){
    print $id;
    
    print $seq;
    print $t;
    print $q;
    $j++;
  }

}
close(INFL);

print STDERR "processed $i reads and got $j mRNA reads\n";

