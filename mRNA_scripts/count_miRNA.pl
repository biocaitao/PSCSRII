use strict;

if(@ARGV<1) {
  die "perl count_miRNA.pl map.anno.deg";
}

#count  all;
my %c1 =();
#multireads
my %mr =();


open(INFL, $ARGV[0]) or die "$!";
while(<INFL>){
#  next if (!/pre/);
  chomp;
  my @t =split /\t/;
  my @t2 = split /_/, $t[0];
  $c1{$t[-1]}{$t2[1]}++;
}
close(INFL);


print "gene\tcell\tcount\n";
foreach my $m (keys %c1){
   foreach my $c (keys %{$c1{$m}}){
     print $m, "\t";
     print $c, "\t";
     print $c1{$m}{$c};
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
