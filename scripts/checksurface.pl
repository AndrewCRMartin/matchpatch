#!/usr/bin/perl -s

use strict;

my $pdbFile = shift(@ARGV);

my $chain = defined($::chain)?$::chain:'A';

my $tmpDir  = "/var/tmp/checksurface_" . $$ . time();
`mkdir $tmpDir`;
die "Can't create $tmpDir directory" if(! -d $tmpDir);

my @surfaceResidues = GetSurfaceResidues($tmpDir, $pdbFile);

foreach my $surfaceResidue (@surfaceResidues)
{
   if($surfaceResidue =~ /^$chain/)
   {
       `pdbmakepatch -r 15 $surfaceResidue CA $tmpDir/access.pdb >$tmpDir/patch.pdb`;
       ExtractPatch($tmpDir, "patch.pdb", "patchonly.pdb");
       `matchpatchsurface -s $tmpDir/patchonly.pdb > $tmpDir/patch_${surfaceResidue}.surf`;
       `pdbgetchain A $pdbFile | matchpatchsurface -s > $tmpDir/patch_pdbfile.surf`;
       `matchpatch $tmpDir/patch_${surfaceResidue}.surf $tmpDir/patch_pdbfile.surf >$tmpDir/results.txt`;
       my $lines = `wc -l $tmpDir/results.txt`;
       if($lines > 2)
       {
           print "\n$surfaceResidue\n";
           system("cat $tmpDir/results.txt");
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
    my($tmpDir, $pdbFile) = @_;
    
    `pdbhetstrip $pdbFile | pdbsolv -x -r $tmpDir/resSolv.dat > $tmpDir/access.pdb`;
    my @surfaceResidues = ();
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
                    push @surfaceResidues, $resid;
                }
            }
        }
    }

    return(@surfaceResidues);
}




