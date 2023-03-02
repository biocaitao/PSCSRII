use strict;

if(@ARGV<1) {
  die "perl count_ncRNA.pl ncRNA.map";
}

#count  all;
my %c1 =();
#multireads
my %mr =();
my %mr2 =();


#gene unique 
my %count1 =();
open(INFL, $ARGV[0]) or die "$!";
while(<INFL>){
  chomp;
  my @t =split /\t/;
  my ($nh) = $t[11] =~  /NH:i:(\d+)/;

  my @t2 = split /_/, $t[0];
  my @t3 = split / /, $t[-1];
  #feach each cell, for each gene;
  foreach my $g (@t3 ){
    $count1{$t2[1]}{$g} += 1/$nh;
  }

}
close(INFL);
#readid gene map 

print "cell\tgene\tcount\n";
foreach my $m (keys %count1){
   foreach my $c (keys %{$count1{$m}}){
     print $m, "\t";
     print $c, "\t";
     print $count1{$m}{$c};
     print "\n"; 
   }
}


sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

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
