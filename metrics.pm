#!/usr/bin/perl -Wall
$| = 1;

use File::Basename;
my $SCRIPTS_DIR = dirname(__FILE__);

my %cmds = ();
#$cmds{"score"}       =  "$SCRIPTS_DIR/score.static.linuxgccrelease"; # -s foo.pdb -ignore_unrecognized_res -corrections::score::score12prime true -corrections::score::no_his_his_pairE
#$cmds{"score"}       =  "/dnas/apps/rosetta/current/source/bin/score.default.linuxgccrelease"; # -s foo.pdb -ignore_unrecognized_res -corrections::score::score12prime true -corrections::score::no_his_his_pairE
#$cmds{"score"}       =  "/dnas/home/dwkulp/software/rosetta_git_073016/Rosetta/main/source/bin/score.linuxgccrelease"; # -s foo.pdb -ignore_unrecognized_res -corrections::score::score12prime true -corrections::score::no_his_his_pairE
$cmds{"score"}       =  "/home/dwkulp/wistar/software/Rosetta/main/source/bin/score.linuxgccrelease";
#$cmds{"rosettadb"}   =  "$SCRIPTS_DIR/rosetta_database";
#$cmds{"rosettadb"}   =  "$SCRIPTS_DIR/rosetta_database_2";
#$cmds{"rosettadb"}   =  "/dnas/apps/rosetta/current/database";
#$cmds{"rosettadb"}   =  "/dnas/home/dwkulp/software/rosetta_git_073016/Rosetta/main/database";
$cmds{"rosettadb"}   =  "/home/dwkulp/wistar/software/Rosetta/main/database";
$cmds{"seqalign"}    =  "$SCRIPTS_DIR/muscle" ; # -i seq.fasta -o align.fasta
$cmds{"cluster"}     =  "$SCRIPTS_DIR/blastclust" ; # -i diff.fasta -p T -L 1 -b T -S 100 -o foo
$cmds{"fasta"}       =  "$SCRIPTS_DIR/printSequence.pl"; # list.txt
$cmds{"unsats"}      =  "$SCRIPTS_DIR/report_buried_unsats_for_plugin.static.linuxgccrelease";
#$cmds{"holes"}       =  "$SCRIPTS_DIR/holes.default.linuxgccrelease";
#$cmds{"holes"}       =  "/dnas/apps/rosetta/current/source/bin/holes.default.linuxgccrelease"; #/dnas/apps/rosetta/current/rosetta_source/bin/holes.linuxgccrelease";
#$cmds{"holes"}       =  "/dnas/home/dwkulp/software/rosetta_git_073016/Rosetta/main/source/bin/holes.linuxgccrelease"; #/dnas/apps/rosetta/current/rosetta_source/bin/holes.linuxgccrelease";
#$cmds{"holes"}       =  "/home/dwkulp/wistar/software/Rosetta/main/source/bin/holes.linuxgccrelease";
$cmds{"holes"}       =  "/home/dwkulp/wistar/software/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease  -parser::protocol holes.xml -holes:dalphaball ~/wistar/software/Rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc"; 
#$cmds{"holes"}       =  "/dnas/home/dwkulp/software/rosetta_git_073016/Rosetta/main/source/bin/holes.linuxgccrelease"; #/dnas/apps/rosetta/current/rosetta_source/bin/holes.linuxgccrelease";
#$cmds{"dalphaball"}  =  "$SCRIPTS_DIR/DAlphaBall.gcc";
#$cmds{"dalphaball"}  =  "/dnas/apps/rosetta/current/source/external/DAlpahBall/DAlphaBall.gcc";
#$cmds{"dalphaball"}  =  "/dnas/apps/rosetta/current/source/external/DAlpahBall/DAlphaBall.gcc";
#$cmds{"dalphaball"}  =  "/home/dwulp/wistar/software/Rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc";
$cmds{"resurfaceSB"} =  "$SCRIPTS_DIR/resurfaceSaltBridges";

# rosetta_scripts.linuxgccrelease  -parser::protocol holes.xml -s test.pdb -holes:dalphaball ~/wistar/software/Rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc 
#/dnas/apps/rosetta/current/source/bin/holes.default.linuxgccrelease  -no_his_his_pairE -overwrite -mute core.conformation -out:file:silent_struct_type binary -database /dnas/apps/rosetta/current/database -s dr-out-H_PIKAAdesTryp-designRelax_run_100.pdb -holes:dalphaball /dnas/apps/rosetta/current/source/external/DAlpahBall/DAlphaBall.gcc

use lib("$SCRIPTS_DIR");
require("pdbline.pm");

@EXPORT = qw"reduceClusters getClusters test getScoreTerm getResfileScoreTerm getHoles getUnsatisfiedPolars";

sub reduceClusters(){
  my $data_hash_ref    = shift;
  my $cluster_by       = shift;
  my $cluster_hash_ref = shift;

  foreach $c (sort keys %{$cluster_hash_ref}){
    my @files = split(",",$$cluster_hash_ref{$c});
    
    # Find the memeber of this cluster with the minimum value
    my $minValue = 1000000;
    my $minMember = "";
    foreach (@files){
      if ($$data_hash_ref{$_}{$cluster_by} < $minValue){
	$minValue = $$data_hash_ref{$_}{$cluster_by};
	$minMember = $_;
      }
    }
    # print "Representative member is $minMember,$minValue\n";

    # Delete the data from all memembers except the representative..
    foreach (@files){
      next if ($_ eq $minMember);

      delete $$data_hash_ref{$_};
    }
  }
  
}

sub getClusters(){
  my $data_hash_ref    = shift;
  my $list             = shift;
  my $cluster_hash_ref = shift;

  # Get array of pdb_names
  open(TMP,$list) or die "Couldn't open $list\n";
  my @pdb_names = <TMP>;
  close(TMP);

  # Create fasta un-aligned
  system("rm tmp.fasta");
  foreach (@pdb_names){
    my $pdb = $_;
    $pdb =~ s/\n//;

    my $cmd = $cmds{"fasta"}." $pdb >> tmp.fasta";
    print "\tCMD: $cmd\n";
    system($cmd);
  }

  # Create multiple-sequence alignment
  $cmd = $cmds{"seqalign"}." -in tmp.fasta -out tmp.aln.fasta > /dev/null 2> /dev/null";
  print "\tCMD: $cmd\n";
  system($cmd);

  # Create clusters
  # Using a similarity threshold of 100 creates a non-redundant set of designs (i.e. doesn't do much clustering)
  # Since there are 671 positions total, two sequences that differ by only 2 residues have a similarity of 669/671= 99.701%. 
  # So, if we wanted to group designs that only differ by 2 residues, we should set our threshold to 99.7. If two sequences
  # differ by 3 residues, their similarity is only 99.55 and they would get grouped to separate clusters.
  # a better way to do this would be to make similarity threshold ( seq. len - 3 ) / seq. length - small constant.
  $cmd = $cmds{"cluster"}." -i tmp.aln.fasta -p T -L 1 -b T -S 100 -o tmp.clusters > /dev/null";
  print "\tCMD: $cmd\n";
  system($cmd);
  
  open(TMP, "tmp.clusters") or die "Couldn't open tmp.clusters: $?\n";
  my $clusterId = 0;
  while ($line = <TMP>){

     my $clusterStr = sprintf("Cluster-%04d", $clusterId++);

     $line =~ s/\n//;
     my @files = split(/\s+/,$line);
     foreach (@files){
       $$data_hash_ref{$_}{"cluster"} = sprintf("%04d",$clusterId);
       $$data_hash_ref{$_}{"members"} = sprintf("%04d",$#files+1);
       $$cluster_hash_ref{$clusterStr} .= $_.","; 
     }
  }

}

# Test for adding data to a hash by-reference
sub test(){
  my $data_hash_ref = shift;
  $$data_hash_ref{"test"} = "foo";
}

sub getScoreTerm(){

  my $data_hash_ref = shift;
  my $pdbs          = shift;
  my $term_list     = shift;
  my $score_file    = shift;

  print "\tReading in score file: $score_file\n";
  # Split the term array..
  my @terms = split(",",$term_list);

  # Create array of pdb file names
  print "Reading in PDBS\n";
  open(TMP, "$pdbs") or die "Couldn't open $pdbs: $?\n";
  my @pdb_names = <TMP>;
  close(TMP);

  # Run scoring app if needed..
  my $ran_score_app=0;
  if (!$score_file){

    system("rm default.sc");
#    my $cmd = "$cmds{\"score\"} -database $cmds{\"rosettadb\"} -l $pdbs -ignore_unrecognized_res -score:weights score12prime.wts -no_his_his_pairE -score::dun10 > /dev/null" ;
    my $cmd = "$cmds{\"score\"}  -database $cmds{\"rosettadb\"} -l $pdbs -ignore_unrecognized_res -no_his_his_pairE > /dev/null" ;
    print "\tCMD: $cmd\n";
    system($cmd);
      # output:
      #SCORE:     score     fa_atr     fa_rep     fa_sol    fa_intra_rep    pro_close    fa_pair    hbond_sr_bb    hbond_lr_bb    hbond_bb_sc    hbond_sc    dslf_ss_dst    dslf_cs_ang    dslf_ss_dih    dslf_ca_dih       rama      omega     fa_dun    p_aa_pp        ref    allatom_rms    gdtmm    gdtmm1_1    gdtmm2_2    gdtmm3_3    gdtmm4_3    gdtmm7_4    irms    maxsub    maxsub2.0    rms       description
      #SCORE:  2547.191  -2354.512   2578.751   1334.732           8.746       65.900    -50.074        -90.993        -66.968        -12.671     -19.297        -18.311          2.512          3.624         10.336    108.569      2.227   1158.601      5.340   -119.320          0.000    1.000       1.000       1.000       1.000       1.000       1.000   0.000   673.000      673.000  0.000 JRFL_chA_ats_0001    

    $score_file = "default.sc";
    $ran_score_app=1;
  }   


  # Parse score file
  open(TMP, "$score_file") or die "Couldn't open $score_file : $?\n";
  my @lines = <TMP>;
  close(TMP);

  # Auto-discover header line, don't assume its line 0 or line 1...
  my @headerLine = grep { $lines[$_] =~ /description/ } 0..$#lines;
  if (@headerLine == 0){
    die "ERROR couldn't find header line defined by 'description'";
  }
  if (@headerLine != 1){
    die "ERROR found more than 1 header line defined by 'description' lines = [@headerLine]";
  }
  #print "HeaderLine is : $lines[$headerLine]\n";


  # Split header into approriate tokens..
  my @header = split(/\s+/,$lines[$headerLine]);

  # Create a local map for each line of header_variables to values
  my %map = ();
  print "Create a map\n";  
  for ($l = $headerLine+1; $l <= $#lines;$l++){
    my @values = split(/\s+/,$lines[$l]);
    for ($i = 0; $i <= $#header; $i++){
      $map{$l}{$header[$i]} = $values[$i];
    }
  }

  # Open a new score file that has the pdb_names inserted in it
  if ($ran_score_app){
      open(NEWSC, ">metrics.sc");
      print NEWSC $lines[$headerLine];
  }

  # Now fill our data_hash_ref with this data..
  print "Fill hash\n";
  my $pdb_names_index = 0;
  for ($l = $headerLine+1; $l <= $#lines;$l++){

    # Use description from score file by default.
    my $description = $map{$l}{description};
    if (exists $map{$l}{path}){
	$description = $map{$l}{path};
	if (substr($map{$l}{path},length($map{$l}{path})-1) ne "/"){
	    $description .= "/";
	}

	$description .= $map{$l}{description};
    }
    
    if ($ran_score_app){
      $description = $pdb_names[$pdb_names_index];
      $description =~ s/\n//;
    }
    $pdb_names_index++;

    # Search for description in pdb_names array..
    my @fileName = grep { $pdb_names[$_] =~ $description } 0..$#pdb_names;
    if (@fileName == 0){
      print "WARNING **** description: $description is not found in list of pdbs\n";
      next;
    }
    if (@fileName != 1){
      print "WARNING *** description: $description is not unique in list of pdbs. matches = [@fileName]\n";
      next;
    }
    $pdb_names[$fileName[0]] =~ s/\n//; # remove end-of-lines

    for ($t = 0; $t <= $#terms;$t++){
      #print "Adding $pdb_names[$fileName[0]],$header[$t],$terms[$t],$map{$l}{$terms[$t]}\n";
      $$data_hash_ref{$pdb_names[$fileName[0]]}{$terms[$t]} = $map{$l}{$terms[$t]};
    }

    if ($ran_score_app){
	for ($i = 0; $i <= $#header; $i++){
	    if ($header[$i] eq "description"){
		print NEWSC $pdb_names[$fileName[0]]." ";
	    } else {
		print NEWSC $map{$l}{$header[$i]}." ";
	    }
        }
        print NEWSC "\n";
    }
			   
  }
  if ($ran_score_app){
    close(NEWSC);
  }

  print "\t\tdone reading score file\n";

}


sub parseResfile() {

  my $resfile       = shift;

  print "\tPARSE RESFILE $resfile\n";
  # Parse resfile
  my %resfile_positions = ();
  open(TMP,"$resfile") or die "Couldn't open file $resfile\n";
  while ($line = <TMP>){
    $line =~ s/\n//;
    if ($line =~ /^start/) {$past_start = 1; next;}
    next if (!$past_start) ;

    my @toks = split(/\s+/,$line);
    $toks[0] =~ s/\s+//g;
    $toks[1] =~ s/\s+//g;
    $resfile_positions{$toks[1]}{$toks[0]} = $toks[2];
    #print "RESFILE: $toks[1],$toks[0],$toks[2].\n";

  }
  close(TMP);

  print "\t\tdone parsing resfile\n";
  return %resfile_positions;
}

sub getResfileScoreTerm(){

  my $data_hash_ref = shift;
  my $pdbs          = shift;
  my $term_list     = shift;
  my $resfile       = shift;

  my %resfile_positions = &parseResfile( $resfile );

  print "\tGetting resfile residue scores...\n";
  # Split the term array..
  my @terms = split(",",$term_list);

  # Create array of pdb file names
  open(TMP, "$pdbs") or die "Couldn't open $pdbs: $?\n";
  my @pdb_names = <TMP>;
  close(TMP);

  my $count = 0;
  my $no_total = $#pdb_names + 1;
  foreach (@pdb_names){
    my $file = $_;
    $count++;
    print "\t\tProcessing file $count of $no_total\n";
    open(TMP, "$file") or die "Couldn't open file $file\n";

    my %header_map = ();
    my %pdb_hash=();
    my $in_score_section = 0;
    my $res_index = 1;
    while ($line = <TMP>){
      $line =~ s/\n//;

      if ($line =~ /^ATOM/ || $line =~ /^HETATM/){

	my $resn = PDBline($line,"resname"); $resn =~ s/\s//g;
	my $resi = PDBline($line,"resseq");$resi =~ s/\s//g;
	my $resc = PDBline($line,"icode");$resc =~ s/\s//g;
	my $chain= PDBline($line,"chainid");$chain =~ s/\s//g;
	my $name = PDBline($line,"atomname");	$name =~ s/\s//g;

	next if ($name ne "CA");

	$pdb_hash{$res_index}{"pdb_chain"} = $chain;
	$pdb_hash{$res_index}{"pdb_resi"} = $resi;
	$pdb_hash{$res_index}{"pdb_resc"} = $resc;
	$pdb_hash{$res_index}{"pdb_resn"} = $resn;

	$res_index++;
	next;
      }

      if ($line =~ /^label/){
	my @header = split(/\s+/,$line);
	# Create a local map for each line of header_variables to values

	for ($i = 0; $i <= $#header; $i++){
	    $header_map{$header[$i]} = $i;
        }

	$in_score_section=1;
	next;
      }

      next if ($line =~ /^weights/);
      next if ($line =~ /^pose/);


      if ($in_score_section && !($line =~ /END_POSE_ENERGIES_TABLE/)){
	my @toks = split(/\s+/,$line);
	my $index = 0;
	if ($toks[0] =~ /_(\d+)$/){
	  $index = $1;
	} else {
	  print "ERROR score line does not have an index: Line: $line\n";
	  exit(234);
	}
	
	my $aa = substr($toks[0],0,3);
	if ($aa ne $pdb_hash{$index}{"pdb_resn"}){
	  if (!($aa eq "CYD" && $pdb_hash{$index}{"pdb_resn"} eq "CYS")){
	    print "ERROR AA($aa) is not the same as ATOM line AA(".$pdb_hash{$index}{"pdb_resn"}.") [ ".$pdb_hash{$index}{"pdb_chain"}.",".$pdb_hash{$index}{"pdb_resi"}."] \n";
	    exit(234);
	  }
	}

	#print "Looking up [ ".$pdb_hash{$index}{"pdb_chain"}.",".$pdb_hash{$index}{"pdb_resi"}."] as $index, $aa value is: ".$resfile_positions{$pdb_hash{$index}{"pdb_chain"}}{$pdb_hash{$index}{"pdb_resi"}}." \n";
	if (! ($resfile_positions{$pdb_hash{$index}{"pdb_chain"}}{$pdb_hash{$index}{"pdb_resi"}} =~ /ALLAA/ ||
	       $resfile_positions{$pdb_hash{$index}{"pdb_chain"}}{$pdb_hash{$index}{"pdb_resi"}} eq "PIKAA" ||
	       $resfile_positions{$pdb_hash{$index}{"pdb_chain"}}{$pdb_hash{$index}{"pdb_resi"}} eq "NATAA") ) {
	  #print "Position was not designable, skip\n";
	  
	  next;
	} else {
	  # print "Position was designable ,add score terms ".$pdb_hash{$index}{"pdb_chain"}." ".$pdb_hash{$index}{"pdb_resi"}."\n";
	}

	
	for ($t = 0; $t <= $#terms;$t++){
	  $toks[$header_map{$terms[$t]}] =~ s/\n//;
	  $file =~ s/\n//;
	  $$data_hash_ref{$file}{"res_".$terms[$t]} += $toks[$header_map{$terms[$t]}];
	}
	
      }

    }
    close(TMP);
  }
  print "\t\tdone getting resfile residue scores...\n";

}


sub getResfileUnsatisfiedPolars() {

    my $data_hash_ref = shift;
    my $pdbs          = shift;
    my $output_file   = shift;
    my $resfile       = shift;

    my %resfile_positions = &parseResfile( $resfile );

    print "\tGetting resfile residue unsatisfied polars...\n";
    # Create array of pdb file names
    open(TMP, "$pdbs") or die "Couldn't open $pdbs: $?\n";
    my @pdb_names = <TMP>;
    close(TMP);

    my $rv;
    my @output;
    if ( !$output_file ) {
        # generate the unsat information on all structures in the list file
#        my $cmd = "$cmds{\"unsats\"} -database $cmds{\"rosettadb\"} -l $pdbs -ignore_unrecognized_res -score::weights score12prime.wts -no_his_his_pairE -overwrite" ;
        my $cmd = "$cmds{\"unsats\"}  -l $pdbs -ignore_unrecognized_res -score::weights score12prime -no_his_his_pairE -overwrite" ;
        print "Running command to get buried unsatisfied polar output...\n$cmd\n";
        @output = `$cmd`;  # doesn't check the return value, but don't know how to catch output and return value with system calls

        $rv = open( UNSATS_OUTPUT, ">metrics.unsats.txt" );
        if ( !defined($rv)) {
            print STDERR "Unable to open file for unsat run output: $!";
            exit;
        }
        foreach my $line ( @output ) {
            print UNSATS_OUTPUT "$line";
        }
        close UNSATS_OUTPUT;

    } else {
        $rv = open( UNSATS_OUTPUT, "$output_file");
        if (!defined($rv)) {
           print STDERR "Unable to open unsat run output file: $!";
           exit;
        }
        @output = <UNSATS_OUTPUT>;
        chomp @output;
        close UNSATS_OUTPUT;
    }

    my %unsat_data;
    foreach my $line ( @output ) {
        # skip lines that don't contain STRUCTURE
        next if ( $line !~ /STRUCTURE/ );
        #print "$line\n";

        # output lines look as follows:
        # apps.pilot.ronj.report_buried_unsats_for_plugin: STRUCTURE:              JRFL_chA_0047   CHAIN:    A   RESIDUE:   711   ATOM:     OXT   SASA:     1.251   HBONDS:    0   BUR_UNSAT:  0

        # set up a regex to find the values
        my $structure;
        my $chain;
        my $residue;
        my $no_buried_unsats;
        if ( $line =~ m/STRUCTURE:\s*([A-Za-z0-9\/\._-]*)\s*CHAIN:\s*([A-Z ])\s*RESIDUE:\s*(\d*).*BUR_UNSAT:\s*([0-9]+)/ ) {
            $structure = $1;
            $chain = $2;
            $residue = $3;
            $no_buried_unsats = $4;
        } else {
            print "Failed to get count of buried polars for residue from output. line: $line\n";
            next;
        }

        #print "Setting hash{ $structure }{ $chain }{ $residue } = $no_buried_unsats\n";
        $unsat_data{ $structure }{ $chain }{ $residue } = $no_buried_unsats;
    }

    # now that we've read the unsats output file into a hash structure, go through all the PDBs and set the number of unsats for just the designed positions in the data hash.
    foreach my $pdb ( @pdb_names ) {
        chomp $pdb;

        foreach my $chain ( keys %resfile_positions ) {
            foreach my $res ( keys %{ $resfile_positions{$chain} } ) {
                #print "resfile_positions{ $chain }{ $res } = $resfile_positions{ $chain }{ $res }\n"; 
                if (! ($resfile_positions{$chain}{$res} eq "ALLAA" ||  $resfile_positions{$chain}{$res} eq "PIKAA" || $resfile_positions{$chain}{$res} eq "ALLAAwc") ){
                    #print "Position was not designable, skip\n";
                    next;
                } else {
                    #print "unsat_data{ $pdb }{ $chain}{ $res } = $unsat_data{ $pdb }{ $chain}{ $res }\n";
                    $$data_hash_ref{ $pdb }{"res_unsats"} += $unsat_data{ $pdb }{ $chain }{ $res };
                }
            }
        }

    }
    print "\t\tdone getting resfile residue unsatisfied polars...\n";

}

sub getUnsatisfiedPolars() {

    my $data_hash_ref = shift;
    my $pdb_list_file = shift;
    my $output_file = shift;

    my $rv;
    my @output;
    if ( !$output_file ) {
        # generate the unsat information on all structures in the list file
#        my $cmd = "$cmds{\"unsats\"} -database $cmds{\"rosettadb\"} -l $pdb_list_file -ignore_unrecognized_res -score::weights score12prime.wts -no_his_his_pairE -overwrite" ;
        my $cmd = "$cmds{\"unsats\"}  -l $pdb_list_file -ignore_unrecognized_res -score::weights score12prime.wts -no_his_his_pairE -overwrite" ;
        print "Running command to get buried unsatisfied polar output...\n$cmd\n";
        @output = `$cmd`;  # doesn't check the return value, but don't know how to catch output and return value with system calls

        $rv = open( UNSATS_OUTPUT, ">metrics.unsats.txt" );
        if ( !defined($rv)) {
            print STDERR "Unable to open file for unsat run output: $!";
            exit;
        }
        foreach my $line ( @output ) {
            print UNSATS_OUTPUT "$line";
        }
        close UNSATS_OUTPUT;

    } else {
        $rv = open( UNSATS_OUTPUT, "$output_file");
        if (!defined($rv)) {
           print STDERR "Unable to open unsat run output file: $!";
           exit;
        }
        @output = <UNSATS_OUTPUT>;
        chomp @output;
        close UNSATS_OUTPUT;
    }

    foreach my $line ( @output ) {
        # skip lines that don't contain STRUCTURE OR INPUT
        #next if ( $line !~ /STRUCTURE/ ) && ( $line !~ /INPUT/ );
        next if ( $line !~ /INPUT/ );
        #print "$line\n";

        # output lines look as follows:
        # apps.pilot.ronj.report_buried_unsats_for_plugin: INPUT:  fixbb_minimize.JRFL_chA.0179_score2   total_polar: 1932   buried_polar: 1010   buried_unsat: 374   exposed_polar: 922   exposed_sat: 127   probe_radius:  1.40
        # for now, all we want is the number of buried unsatisfied polars for each structure
        # later, we'll want to get the per-residue counts
        #my @fields = split( /\s+/, $line );
        #my $structure = $fields[2];  # note, this doesn't include the .pdb extension
        #my $no_buried_unsats = $fields[8];

        # can't do it by splitting the line on whitespace because the PDB name and counts are in different columns
        # depending on if the executable has MPI-enabled
        # instead set up a regex to find the values
        my $structure;
        my $no_buried_unsats;
        if ( $line =~ m/INPUT:\s*([A-Za-z0-9\/\._-]*)\s*total_polar:\s*\d*\s*buried_polar:\s*(\d*)/ ) {
            $structure = $1;
            $no_buried_unsats = $2;
        } else {
            print "Failed to get count of buried polars from output. line: $line\n";
            next;
        }
 
        #print "Setting hash{ $structure }->{ no. buried unsats } = $no_buried_unsats\n";
        $$data_hash_ref{ $structure }{ "unsats" } = $no_buried_unsats
    }

}
sub getSaltBridgeScores() {
    my $data_hash_ref = shift;
    my $pdb_list_file = shift;

    # Get array of pdb_names
    open(TMP,$pdb_list_file) or die "Couldn't open $pdb_list_file\n";
    my @pdb_names = <TMP>;
    close(TMP);

    foreach (@pdb_names){
	my $pdb = $_;
	$pdb =~ s/\n//;
	my $cmd = "$cmds{\"resurfaceSB\"} --pdb $pdb --scoreOnly 2> /dev/null | grep SCORE";
	open(TMP, "$cmd |");
	my $sb_score  = 0;
	while ($line = <TMP>){
	    $line =~ s/\n//g;
	    if ($line =~ /SCORE:\s+(\S+)/){
		$sb_score = $1;
	    }
	}
    
	# Add $pdb $sb_score to $$data_hash_ref
	$$data_hash_ref{$pdb}{"saltBridgeScore"} = $sb_score;

	close(TMP);
    }

    
}

sub getHoles() {

    my $data_hash_ref = shift;
    my $pdb_list_file = shift;
    my $output_file = shift;

    my $rv;
    my @output;
    if ( !$output_file ) {
        # generate the holes information on all structures in the list file
#        my $cmd = "$cmds{\"holes\"} -database $cmds{\"rosettadb\"} -l $pdb_list_file -holes:dalphaball $cmds{\"dalphaball\"} -ignore_unrecognized_res -score::weights score12prime.wts -no_his_his_pairE -overwrite -mute core.conformation -out:file:silent_struct_type binary" ;
        my $cmd = "$cmds{\"holes\"}  -l $pdb_list_file -ignore_unrecognized_res -overwrite -mute core.conformation -out:file:silent_struct_type binary" ;
        print "Running command to get packing output...\n$cmd\n";
        @output = `$cmd`;  # doesn't check the return value, but don't know how to catch output and return value with system calls

        $rv = open( HOLES_OUTPUT, ">metrics.holes.txt" );
        if ( !defined($rv)) {
            print STDERR "Unable to open file for holes run output: $!";
            exit;
        }
        foreach my $line ( @output ) {
            print HOLES_OUTPUT "$line";
        }
        close HOLES_OUTPUT;

    } else {
	print "READING HOLES FILE $output_file\n";
        $rv = open( HOLES_OUTPUT, "$output_file");
        if (!defined($rv)) {
           print STDERR "Unable to open holes run output file: $!";
           exit;
        }
        @output = <HOLES_OUTPUT>;
        chomp @output;
        close HOLES_OUTPUT;
    }

#protocols.jd2.PDBJobInputter: PDBJobInputter::pose_from_job
#protocols.jd2.PDBJobInputter: filling pose from saved copy 5VK2_LASV_trimer_clean_0007.pdb
#protocols.rosetta_scripts.ParsedProtocol: =======================BEGIN FILTER full_holes=======================
#protocols.rosetta_scripts.ParsedProtocol: =======================END FILTER full_holes=======================
#protocols.rosetta_scripts.ParsedProtocol: setting status to success
#protocols.simple_filters.HolesFilter:  computing using full pose 
#core.scoring.packing.compute_holes_score: compute_holes_surfs try: 1
#core.scoring.packing.compute_holes_score: compute_holes_surfs completed successfully
#core.scoring.packing.compute_holes_score: compute_rosettaholes_score done: 1.88266

    foreach my $line ( @output ) {
        chomp $line;
        #print "$line\n";
#        next if ( $line !~ /apps\.RosettaHoles/ );  # only look at the app output lines
        next if ( $line !~ /^RosettaHoles/ );  # only look at the app output lines
        next if ( $line =~ /HEADER/ );  # skip the HEADER line

        # output lines look as follows:
        # apps.RosettaHoles:              JRFL_chA_0001.pdb  5.9703  2.9799  1.6782  5.0874 15852.95   8927.79  27065.00   314.382

        # set up a regex to find the values
        my $structure;
        my $holes_score;
	my $holes_dec;
#        RosettaHoles: HEADER         File/Tag         RosHol    resl   dec25   dec15   AAresl   AAdec25   AAdec15  RosScore
#        RosettaHoles:        4fp8_C05_cut2_renum.pdb  1.6573  1.6553 -1.7870 -1.9599  2136.96  -2307.08  -2530.28  -180.029

#        if ( $line =~ m/RosettaHoles:\s*([A-Za-z0-9\/\._-]*)\s*[0-9]*\.?[0-9]+\s*[0-9]*\.?[0-9]+\s*[0-9]*\.?[0-9]*\s*([0-9]*\.?[0-9]+)/ ) {
        if ( $line =~ m/RosettaHoles:\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+/){
            $structure   = $1;
	    $holes_score = $2;
            $holes_dec   = $3;
        } else {
            print "Failed to get rosetta holes score from output. line: $line\n";
            next;
        }

        #print "Setting hash{ $structure }->{ holes} = $holes_score\n";
        $$data_hash_ref{ $structure }{ "holes" }       = $holes_score;
	$$data_hash_ref{ $structure }{ "holes_decoy" } = $holes_dec;
    }

}

