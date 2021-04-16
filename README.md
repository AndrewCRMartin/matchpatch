MatchPatch
==========

Match patches of proteins using the Lesk graph matching algorithm (see
Brint and Willett (1986) 'Pharmacophore Pattern Matching In Files Of
3D Chemical Structures: Comparison OfGeometric Searching Algorithms',
*J Mol Graph* **5**:49-56. The method is rather slow, but simple to implement.

To use the program you first generate the patterns using the `surface` program. e.g.

```
surface file.pdb > file.surf
```

(optionally a file can be provided that specifies residue ranges to be considered).

This identifies surface residues that are charged or aromatic and
outputs pairs of these residues with the distance between them.

**This must be done to generate the pattern for which you wish to search and for the protein against which you wish to search**.

The `match` program then implements Lesk's algorithm to compare the
pattern with the protein:

```
match pattern.surf protein.surf
```


History
-------

This code was originally written while working at the DKfz, Heidelberg
in 1993. It was updated in 2021 to compile cleanly using newer compilers
and to call the new BiopLib routines.
