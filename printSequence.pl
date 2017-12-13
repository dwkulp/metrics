#!/usr/bin/perl
my @inputs = @ARGV;
if($#inputs < 0) {
   die "Need a pdb file!\n";
}

my %resCode;
$resCode{"ALA"} = "A";
$resCode{"CYS"} = "C";
$resCode{"ASP"} = "D";
$resCode{"GLU"} = "E";
$resCode{"PHE"} = "F";
$resCode{"GLY"} = "G";
$resCode{"HIS"} = "H";
$resCode{"HSD"} = "H";
$resCode{"HSE"} = "H";
$resCode{"ILE"} = "I";
$resCode{"LYS"} = "K";
$resCode{"LEU"} = "L";
$resCode{"MET"} = "M";
$resCode{"ASN"} = "N";
$resCode{"PRO"} = "P";
$resCode{"GLN"} = "Q";
$resCode{"ARG"} = "R";
$resCode{"SER"} = "S";
$resCode{"THR"} = "T";
$resCode{"VAL"} = "V";
$resCode{"TRP"} = "W";
$resCode{"TYR"} = "Y";

unless (open PROTEINS, $inputs[0]) {
   die "Cannot open logfile: $!";
}
my $chain = "";
my $resNum = "999";
my $sequence = "";
my @seqs;
my @chains;
while (<PROTEINS>) {
   chomp;
   my $entry = $_;
   if (substr($entry,0,4) eq "ATOM") {
      my $currentRes = substr($entry,17,3);
      my $currentChain = substr($entry,21,1);
      my $currentResNum = substr($entry,23,3);
      my $currentSegid  = substr($entry,72,4);

      #print "CURSEGID: $currentSegid.\n";
      next if ($inputs[1] && "$currentSegid" ne "HHHH");

      if ($currentChain ne $chain) {
         if ($chain ne "") { push(@seqs, $sequence); }
         push(@chains, $currentChain);
         $chain = $currentChain;
         $sequence = $resCode{$currentRes};
         $resNum = $currentResNum;
      }
      elsif ($resNum ne $currentResNum) {
         $sequence .= $resCode{$currentRes};
         $resNum = $currentResNum;
      }
   }
}
push(@seqs, $sequence);
close PROTEINS;

print ">$inputs[0]\n";
for (my $i = 0; $i <= $#seqs; $i++) {
   #print $chains[$i] . "\t" . $seqs[$i] . "\n";
  print $seqs[$i] . "\n";
}

