#!/usr/bin/perl
$| = 1;

use FindBin;
use lib $FindBin::Bin;

require "metrics.pm";

use Getopt::Long;
&GetOptions(
	"list=s"                      => \$list,
	"cluster"                     => \$cluster,
	"reduce_cluster_by=s"         => \$cluster_by,
	"scoreterms=s"                => \$scoreterms,
	"scorefile=s"                 => \$scorefile,
	"buried_unsats"               => \$buried_unsats,
	"buried_unsats_output_file=s" => \$buried_unsats_output_file,
	"holes"                       => \$holes,
	"holes_output_file=s"         => \$holes_output_file,
	"resfile=s"                   => \$resfile,
        "sb_score"                    => \$sb_score,
	"debug"                       => \$debug,
	"help"                        => \$help_flag,
);



if (!$list) {
  print "ERROR **** Need to define a list\n";
  $help_flag=1;
}

if (!$scoreterms){
  print "ERROR **** Need to define scoreterms\n";
  $help_flag=1;
}

if ($cluster_by && !($scoreterms =~ /$cluster_by/)){	
  print "WARNING cluster_by($cluster_by) needs to be in scoreterms($scoreterms), adding it automatically\n";
  $scoreterms .= ",$cluster_by";
}

if ($help_flag) {
print <<OUTEND;

Usage:
getMetrics.pl --list LIST_OF_DESIGNED_PDBS
    --cluster
        --reduce_cluster_by score_type
    --scoreterms "fa_atr,rama,fa_dun,ref"
        --scorefile default.sc
    --buried_unsats
        --buried_unsats_output_file output.txt
    --holes
        --holes_output_file output.txt
    --sb_score
    --resfile resfile
    --debug
    --help
         
OUTEND
exit();
}

# Container for the structure-specific data
my %struct_data = ();

# Container for the residue-specific data
my %resi_data   = ();

# Here is a test
if ($debug) {
  test(\%struct_data);
  foreach $k (keys %struct_data){
    print "$k in data is $struct_data{$k}\n";
  }
  die "Done.";
}


if ($scoreterms){
  getScoreTerm(\%struct_data,$list,$scoreterms,$scorefile);

  print "Done: getScoreTerm\n";
  select(STDERR);
  $| = 1;
	  
  select(STDOUT);
  $| = 1;
}

if ($resfile){
  getResfileScoreTerm(\%struct_data,$list,$scoreterms,$resfile);
  print "Done: getResfileScoreTerm\n";
  if ($buried_unsats) {
    getResfileUnsatisfiedPolars(\%struct_data,$list,$buried_unsats_output_file,$resfile);
    print "Done: getResfileUnsatisfiedPolars\n";
  }
  select(STDERR);
  $| = 1;
	  
  select(STDOUT);
  $| = 1;
}

if ($buried_unsats){
  getUnsatisfiedPolars(\%struct_data,$list,$buried_unsats_output_file);
  print "Done: getUnsatisfiedPolars\n";
}

if ($holes){
  getHoles(\%struct_data,$list,$holes_output_file);
  print "Done: getHoles\n";
}

if ($sb_score){
    getSaltBridgeScores(\%struct_data,$list);
    print "Done: getSaltBridgeScores\n";
}

if ($cluster){
  my %clusters = ();
  getClusters(\%struct_data,$list,\%clusters);
  print "Done: getClusters\n";

  # Reduce each cluster to top scoring 
  if ($cluster_by){
    reduceClusters(\%struct_data,$cluster_by,\%clusters);
    print "Done: reduceClusters\n";
  }
  select(STDERR);
  $| = 1;
	  
  select(STDOUT);
  $| = 1;
}


print "\n\n *** SUMMARY *** \n\n";

# Print out header
my $maxFileNameLength = 0;
my $maxFieldLength = 0;
foreach $d (sort keys %struct_data){
  if (length($d)+6 > $maxFileNameLength){
    $maxFileNameLength = length($d)+6;
  }
  foreach $m (sort keys %{$struct_data{$d}}){
    if (length($m)+4 > $maxFieldLength){
      $maxFieldLength = length($m)+4;
    }
  }
}
foreach $d (sort keys %struct_data){

  print sprintf("%-".$maxFileNameLength."s ","FILE ");
  foreach $m (sort keys %{$struct_data{$d}}){
    print sprintf("%".$maxFieldLength."s   ",$m);
  }
  print "\n";
  last;
}

# Print out data
foreach $d (sort keys %struct_data){
  next if ($d eq "clusters") ;

  print sprintf("%-".$maxFileNameLength."s ",$d);
  foreach $m (sort keys %{$struct_data{$d}}){
	  print sprintf("%".$maxFieldLength.".2f   ",$struct_data{$d}{$m});
  }
  print "\n";
}




print "\n\n";


