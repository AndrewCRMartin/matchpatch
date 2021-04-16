MatchPatch
==========

Match patches of proteins using the Lesk graph matching algorithm (see
Brint and Willett (1986) 'Pharmacophore Pattern Matching In Files Of
3D Chemical Structures: Comparison Of Geometric Searching Algorithms',
*J Mol Graph* **5**:49-56). The method is rather slow, but simple to
implement.

To use the program you first generate the patterns using the `surface`
program. e.g.

```
surface file.pdb > file.surf
```

(optionally a file can be provided that specifies residue ranges to be
considered).

This identifies surface residues that are charged or aromatic and
outputs pairs of these residues with the distance between them.

**This must be done to generate the pattern for which you wish to search and for the protein against which you wish to search**.

The `match` program then implements Lesk's algorithm to compare the
pattern with the protein:

```
match pattern.surf protein.surf
```

Type `surface -h` or `match -h` for help.

Compiling
---------

You must first install BiopLib
(https://github.com/ACRMGroup/bioplib/). If this is installed using
the default location, you then simply need to type:

```
make
make install
```

This will install `surface` and `match` in your ~/bin directory
(creating it if it doesn't exist).


History
-------

This code was originally written while working at the DKfz, Heidelberg
in 1993. It was updated in 2021 to compile cleanly using newer compilers
and to call the new BiopLib routines.
