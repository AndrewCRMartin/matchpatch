#!/usr/bin/perl -s

use strict;

my $chain    = defined($::chain)?$::chain:'A';
my $minMatch = defined($::min)?$::min:3;
my $probe    = defined($::probe)?$::probe:1.4;

UsageDie($chain, $minMatch) if((scalar(@ARGV) < 2) || defined($::h));

my $templatePdbFile = shift(@ARGV);

my $tmpDir  = "/var/tmp/checksurface_" . $$ . time();
`mkdir $tmpDir`;
die "Can't create $tmpDir directory" if(! -d $tmpDir);

# Get a list of surface residues from the template PDB file
my @templateSurfaceResidues =
    GetSurfaceResidues($tmpDir, $templatePdbFile, $probe);

foreach my $targetPDBFile (@ARGV)
{
    # Extract the surface descriptor for this PDB file
    `matchpatchsurface $targetPDBFile > $tmpDir/targetPDBFile.surf`;

    foreach my $templateSurfaceResidue (@templateSurfaceResidues)
    {
        # If the residue is in the chain of interest
        if($templateSurfaceResidue =~ /^$chain/)
        {
            # Build a patch around that residue
            `pdbmakepatch -r 15 $templateSurfaceResidue CA $tmpDir/access.pdb >$tmpDir/patch.pdb`;
            ExtractPatch($tmpDir, "patch.pdb", "patchonly.pdb");
            
            # Create the surface descriptor for the patch
            if(! -f "$tmpDir/patch_${templateSurfaceResidue}.surf")
            {
                `matchpatchsurface -s $tmpDir/patchonly.pdb > $tmpDir/patch_${templateSurfaceResidue}.surf`;
            }

            # Match this patch against the surface of the PDB file
            `matchpatch $tmpDir/patch_${templateSurfaceResidue}.surf $tmpDir/targetPDBFile.surf >$tmpDir/results.txt`;


            my $lines = `wc -l $tmpDir/results.txt`;
            if($lines >= $minMatch)
            {
                print "\nTemplate patch around $templateSurfaceResidue matches $targetPDBFile:\n";
                system("cat $tmpDir/results.txt");
            }
        }
    }
}



unlink $tmpDir;



sub ExtractPatch
{
    my($tmpDir, $inPDB, $outPDB) = @_;

    if(open(my $in, '<', "$tmpDir/$inPDB"))
    {
        if(open(my $out, '>', "$tmpDir/$outPDB"))
        {
            while(<$in>)
            {
                if(substr($_, 62, 4) eq "1.00")
                {
                    print $out "$_";
                }
            }
            close($out);
        }
        close($in);
    }
}



sub GetSurfaceResidues
{
    my($tmpDir, $templatePdbFile, $probe) = @_;
    
    `pdbhetstrip $templatePdbFile | pdbsolv -p $probe -x -r $tmpDir/resSolv.dat > $tmpDir/access.pdb`;
    my @templateSurfaceResidues = ();
    if(open(my $fp, '<', "$tmpDir/resSolv.dat"))
    {
        while(<$fp>)
        {
            if(/^RESACC/)
            {
                my(@fields) = split;
                if($fields[7] > 10.0)
                {
                    my $resid = $fields[1] . $fields[2];
                    push @templateSurfaceResidues, $resid;
                }
            }
        }
    }

    return(@templateSurfaceResidues);
}


sub UsageDie
{
    my($chain, $minMatch) = @_;
    
    print <<__EOF;

checksurface V1.0 (c) 2021 UCL, Prof Andrew C.R. Martin

Usage: checksurface [-chain=c][-min=n] template.pdb target.pdb [target.pdb ...]
          -chain  Specify the chain of interest in the template [$chain]
          -min    Minimum number of residues to match in a pattern [$minMatch]

Takes a template PDB file and identifies the surface residues. Each of these is used as
the centre of a 15A radius surface patch. For each of these patches, a matchpatch 
surface descriptor is created scanned against the surface of each target PDB file in 
turn.

Note that `matchpatch`, `matchpatchsurface` and BiopTools must be installed and in your
path.
    
__EOF
    
    exit 0;
}

