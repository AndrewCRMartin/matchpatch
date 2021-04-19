/*************************************************************************

   Program:    match
   File:       match.c
   
   Version:    V2.0
   Date:       16.04.21
   Function:   Match 2 distance matrices as created by surface
   
   Copyright:  (c) SciTech Software / abYinformatics 1993-2021
   Author:     Prof. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without print permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and 
   that the source code is provided with the program.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  21.11.93 Original
   V1.1  16.04.21 Cleaned up for modern compilers and BiopLib2
   V1.2  16.04.21 Made inversion of one pattern an option (instead of
                  always doing it)
   V1.3  16.04.21 Now uses our standard way of parsing the command line
   V2.0  16.04.21 Now reads coordinates and features from input files
                  instead of distance matrix which is calculated here

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/fsscanf.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines
*/
#define MAXPROP         8
#define PROP_POSITIVE   0
#define PROP_NEGATIVE   1
#define PROP_AROM       2
#define MAXITER       100
#define MAXDIST        32
#define MAXBUFF       160

#define DEFBIN        1.0     /* Default distance bin size              */
#define DEFACC       50.0     /* Default string match accuracy          */

/************************************************************************/
/* Structure and type definitions
*/
typedef struct
{
   int  resnum[2];
   char resnam[2][8],
        chain[2][8],
        insert[2][8],
        resid[2][16];
   int  dist;
   BOOL dead;
}  DATA;

typedef struct _indata
{
   REAL x, y, z;
   int charge, aromatic, hydropathy;
   char resnam[8], resid[16];
   struct _indata *next;
}  INDATA;

   

typedef struct
{
   int  resnum;
   char resnam[8],
        chain[8],
      insert[8],
      resid[16];
   char dist[MAXDIST],
        property[MAXPROP];
}  ATOM;

/************************************************************************/
/* Globals
*/
REAL gBin      = DEFBIN,    /* Bin size for distance matrix             */
     gAccuracy = DEFACC;    /* Percentage accuracy for string comparison*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *PatFile, char *StrucFile,
                  char *outfile, BOOL *invert);
void Usage(void);
void MatchFiles(FILE *out, FILE *fp_pat, FILE *fp_struc, BOOL invert);
DATA *ReadMatrix(FILE *fp, int *outnatom);
ATOM *CreateAtomArray(DATA *data, int ndata, int *outnatom,
                      BOOL SwapProp);
int ConvertDistanceToBin(REAL dist);
int GotAtom(ATOM *outdata, int natom, char *resid);
void FillAtom(ATOM *outdata, int natom, int pos, char *resid,
              char *resnam, int DistRange, BOOL SwapProp);
void DoLesk(FILE *out, int npat, DATA *pat, int nstruc, DATA *struc,
            BOOL invert);
void KillAtom(char *resid, DATA *data, 
              int ndata);
void TrimBitStrings(int npat, ATOM *pat, int nstruc, ATOM *struc);
void PrintResults(FILE *out, int NPatAtom,   ATOM *PatAtom, 
                  int NStrucAtom, ATOM *StrucAtom);
void PrintBestMatch(FILE *out,
                    ATOM *PatAtom,   int PatIndex, 
                    ATOM *StrucAtom, int NStrucAtom);
BOOL Compare(char *str1, char *str2, int len);
REAL CalcScore(char *str1, char *str2, int len);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for matching output files from surface.

   18.11.93 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char PatFile[MAXBUFF],
        StrucFile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *fp_pat = NULL,
        *fp_struc = NULL,
        *out = stdout;
   BOOL invert;

   if(ParseCmdLine(argc, argv, PatFile, StrucFile, outfile, &invert))
   {
      if((fp_pat = fopen(PatFile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open pattern file: %s\n",PatFile);
         exit(1);
      }
      if((fp_struc = fopen(StrucFile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open pattern file: %s\n",StrucFile);
         exit(1);
      }
      if(outfile[0] && ((out = fopen(outfile, "w"))==NULL))
      {
         fprintf(stderr,"Unable to write outpuf file: %s\n",outfile);
         exit(1);
      }

      MatchFiles(out, fp_pat, fp_struc, invert);
   }
   else
   {
      Usage();
   }
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *PatFile, 
                     char *StrucFile, char *outfile, BOOL *invert)
   ---------------------------------------------------------------
   Read the command line

   18.11.93 Original   By: ACRM
   22.11.93 Added flags
   16.04.21 Added -i flag
   16.04.21 Rewritten
*/
BOOL ParseCmdLine(int argc, char **argv, char *PatFile, char *StrucFile,
                  char *outfile, BOOL *invert)
{
   argc--;
   argv++;
   
   PatFile[0]  = StrucFile[0] = outfile[0] = '\0';
   *invert = FALSE;
   
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd': case 'D':
            argc--; argv++;
            sscanf(argv[0],"%lf",&gBin);
            if(gBin == 0.0) gBin = 1.0;
            break;
         case 'a': case 'A':
            argc--; argv++;
            sscanf(argv[0],"%lf",&gAccuracy);
            if(gAccuracy == 0.0) gAccuracy = 100.0;
            break;
         case 'i': case 'I':
            *invert = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
         argc--;
         argv++;
      }
      else
      {
         /* Check that there are 2-3 arguments left                     */
         if((argc < 2) || (argc > 3))
            return(FALSE);

         strcpy(PatFile,argv[0]);
         argc--; argv++;
         strcpy(StrucFile,argv[0]);
         argc--; argv++;

         if(argc)
         {
            strcpy(outfile, argv[0]);
            argc--; argv++;
         }
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Display a usage message

   18.11.93 Original   By: ACRM
   22.11.93 Added flag decriptions
   16.04.21 V1.1, V1.2, V1.3, V2.0
*/
void Usage(void)
{
   fprintf(stderr,"\nMatch V2.0 (c) 1993-2021 SciTech Software / \
abYinformatics\n");

   fprintf(stderr,"\nUsage: match [-i][-d binsize][-a accuracy] \
patternFile structureFile [outfile]\n");
   fprintf(stderr,"       -i invert the properties in the pattern \
file\n");
   fprintf(stderr,"       -d specifies distance bin size \
(default: %.1f)\n", (double)DEFBIN);
   fprintf(stderr,"       -a specifies percent dist string match \
accuracy (default: %.1f)\n", (double)DEFACC);
   fprintf(stderr,"\nFind potential matches for a pattern in a structure \
using Lesk's method\n");
   fprintf(stderr,"The input files are generated by `surface`\n");
   fprintf(stderr,"Currently handles properties for charged and aromatic \
residues\n\n");
}

/************************************************************************/
/*>void MatchFiles(FILE *out, FILE *fp_pat, FILE *fp_struc, BOOL invert)
   ---------------------------------------------------------------------
   18.11.93 Original   By: ACRM
*/
void MatchFiles(FILE *out, FILE *fp_pat, FILE *fp_struc, BOOL invert)
{
   DATA *pat,
        *struc;
   int  npat,
        nstruc;

   pat   = ReadMatrix(fp_pat,   &npat);
   struc = ReadMatrix(fp_struc, &nstruc);

   fprintf(stderr, "%d records read from pattern; \
%d records read from structure\n",
           npat, nstruc);

   DoLesk(out, npat, pat, nstruc, struc, invert);
}

/************************************************************************/
/*>ATOM *ReadMatrix(FILE *fp, int *natom)
   --------------------------------------
   Read the output from surface. Create an array of type ATOM which
   contains flags for the presence of connections of a given distance and
   flags the atom type.

   18.11.93 Original   By: ACRM
   22.11.93 Corrected return values
   19.04.21 Now reads the coordinates and does the distance calculations
*/
DATA *ReadMatrix(FILE *fp, int *outnatom)
{
   int  maxrec   = 0,
        i        = 0;
   DATA *outdata = NULL;
   char buffer[MAXBUFF];
   INDATA *indata = NULL,
          *ini    = NULL,
          *inj    = NULL;
   
   *outnatom = 0;

   /* Scan the file to see how long it is                               */
   while(fgets(buffer,MAXBUFF-1,fp)) maxrec++;
   rewind(fp);

   /* maxrec is the number of atom points. We need the number of distances
      between these atoms which is (maxrec^2 - maxrec/2) - one 
      off-diagonal triangle from the matrix
   */
   maxrec = (int)((REAL)(maxrec*maxrec)/2.0 - (REAL)maxrec/2.0);

   /* Allocate this much space                                          */
   if((outdata = malloc(maxrec * sizeof(DATA)))==NULL)
      return(NULL);

   /* Now read the file again, filling in data                          */
   i=0;
   while(fgets(buffer,MAXBUFF-1,fp))
   {
      if(indata==NULL)
      {
         INIT(indata, INDATA);
         ini = indata;
      }
      else
      {
         ALLOCNEXT(ini, INDATA);
      }

      if(ini == NULL)
      {
         fprintf(stderr,"No memory for input data\n");
         FREELIST(indata, INDATA);
         exit(1);
      }
      sscanf(buffer,"%s %s %lf %lf %lf %d %d %d",
             ini->resnam, ini->resid,
             &ini->x, &ini->y, &ini->z,
             &ini->charge, &ini->aromatic, &ini->hydropathy);
   }

   i = 0;
   for(ini=indata; ini!=NULL; NEXT(ini))
   {
      for(inj=ini->next; inj!=NULL; NEXT(inj))
      {

         int  DistRange;
         REAL dist;
         
         dist = DIST(ini, inj);
         DistRange = ConvertDistanceToBin(dist);

         
         strcpy(outdata[i].resnam[0], ini->resnam);
         strcpy(outdata[i].resid[0], ini->resid);

         strcpy(outdata[i].resnam[1], inj->resnam);
         strcpy(outdata[i].resid[1], inj->resid);

         outdata[i].dist = DistRange;
         outdata[i].dead = FALSE;

         i++;
      }
      
   }

   *outnatom = i;
   return(outdata);
}

/************************************************************************/
/*>ATOM *CreateAtomArray(DATA *data, int ndata, int *outnatom, 
                         BOOL SwapProp)
   -----------------------------------------------------------
   19.11.93 Original   By: ACRM
*/
ATOM *CreateAtomArray(DATA *data, int ndata, int *outnatom, BOOL SwapProp)
{
   int  i,
        natom = 0,
        pos;
   ATOM *outatom;

   /* Allocate memory for the output atom array                         */
   if((outatom = (ATOM *)malloc(ndata * sizeof(ATOM))) == NULL)
      return(NULL);

   for(i=0; i<ndata; i++)
   {
      /* Skip this record if the dead flag is set                       */
      if(data[i].dead) continue;

      /* See if we've got a record for the first residue                */
      pos = GotAtom(outatom, natom, data[i].resid[0]);
      if(pos == (-1)) pos = natom++;

      /* Fill this into the data array                                  */
      FillAtom(outatom, natom, pos,
               data[i].resid[0],
               data[i].resnam[0], data[i].dist,
               SwapProp);

      /* See if we've got a record for the second residue               */
      pos = GotAtom(outatom, natom, data[i].resid[1]);
      if(pos == (-1)) pos = natom++;

      /* Fill this into the data array                                  */
      FillAtom(outatom, natom, pos,
               data[i].resid[1],
               data[i].resnam[1], data[i].dist,
               SwapProp);
   }

   *outnatom = natom;

   return(outatom);
}

/************************************************************************/
/*>int ConvertDistanceToBin(REAL dist)
   -----------------------
   Converts a distance (REAL) to an integer bin number < MAXDIST. The bin
   size is read from the global variable gBin

   19.11.93 Original   By: ACRM
   22.11.93 Changed to read bin size from global gBin
*/
int ConvertDistanceToBin(REAL dist)
{
   int idist;

   idist = (int)(dist/gBin);
   if(idist >= MAXDIST-1) idist = MAXDIST-2;

   return(idist);
}

/************************************************************************/
/*>int GotAtom(ATOM *outdata, int natom, char *resid)
   --------------------------------------------------------------
   Searches the current outdata array to see if we already have this
   residue. If so, returns the index into the array. If not, returns -1

   19.11.93 Original   By: ACRM
   19.04.21 Changed to use resid
*/
int GotAtom(ATOM *outdata, int natom, char *resid)
{
   int i;

   for(i=0; i<natom; i++)
   {
      if(!strcmp(outdata[i].resid, resid))
      {
         return(i);
      }
   }
   return(-1);
}

/************************************************************************/
/*>void FillAtom(ATOM *outdata, int natom, int pos, char *resid,
                 char *resnam, int DistRange,
                 BOOL SwapProp)
   --------------------------------------------------------------------
   Fill in an item in the data array. If (pos == natom-1) then it's a
   new residue so we must fill in all data; otherwise just set the
   appropriate flags.

   19.11.93 Original   By: ACRM
   22.11.93 Added aromatic support
   19.05.94 Added DNA support
   19.04.21 Changed to resid
*/
void FillAtom(ATOM *outdata, int natom, int pos,  char *resid,
              char *resnam, int DistRange, BOOL SwapProp)
{
   int i;

   if(pos == natom-1)
   {
      /* A new residue; fill in all data                                */
      strcpy(outdata[pos].resid,  resid);
      strcpy(outdata[pos].resnam, resnam);

      /* Clear all distance and property flags                          */
      for(i=0; i<MAXPROP; i++) outdata[pos].property[i] = '0';
      outdata[pos].property[MAXPROP-1] = '\0';
      for(i=0; i<MAXDIST; i++) outdata[pos].dist[i]     = '0';
      outdata[pos].dist[MAXDIST-1] = '\0';

      /* Set the property flags for this residue                        */
      if(!strncmp(resnam,"LYS",3) || !strncmp(resnam,"ARG",3))
         outdata[pos].property[SwapProp?PROP_NEGATIVE:PROP_POSITIVE]
            = '1';
      if(!strncmp(resnam,"GLU",3) || 
         !strncmp(resnam,"ASP",3) ||
         !strncmp(resnam,"  A",3) ||
         !strncmp(resnam,"  T",3) ||
         !strncmp(resnam,"  C",3) ||
         !strncmp(resnam,"  G",3))
         outdata[pos].property[SwapProp?PROP_POSITIVE:PROP_NEGATIVE]
            = '1';
      if(!strncmp(resnam,"TYR",3) || 
         !strncmp(resnam,"PHE",3) ||
         !strncmp(resnam,"TRP",3))
         outdata[pos].property[PROP_AROM] = '1';
   }

   /* Now set the distance flag                                         */
   outdata[pos].dist[DistRange] = '1';
}

/************************************************************************/
/*>void DoLesk(FILE *out, int npat, DATA *pat, int nstruc, DATA *struc)
   ---------------------------------------------------------
   Does the actual Lesk pattern matching algorithm (with some 
   modifications).

   19.11.93 Original   By: ACRM
   21.11.93 Added property comparison and printing of results :-)
   16.04.21 Added invert parameter instead of always inverting the pattern
*/
void DoLesk(FILE *out, int npat, DATA *pat, int nstruc, DATA *struc,
            BOOL invert)
{
   ATOM *PatAtom,
        *StrucAtom;
   int  NPatAtom,
        NStrucAtom,
        i, j, k,
        PrevPatAtoms   = 0,
        PrevStrucAtoms = 0;
   
   for(i=0; i<MAXITER; i++)
   {
      /* Create the bit strings for the atoms from the data arrays      */
      PatAtom   = CreateAtomArray(pat,   npat,   &NPatAtom,   invert);
      StrucAtom = CreateAtomArray(struc, nstruc, &NStrucAtom, FALSE);

      /* Print information on remaining atoms                           */
      fprintf(stderr, "Iteration %d: %d pattern atoms and %d structure \
atoms remain\n",i,NPatAtom,NStrucAtom);

      /* Remove any distance flags from the bit strings which are 
         never seen in the other structure
      */
      TrimBitStrings(NPatAtom, PatAtom, NStrucAtom, StrucAtom);

      /* Exit if we've converged                                        */
      if(NPatAtom == PrevPatAtoms && NStrucAtom == PrevStrucAtoms)
         break;
      PrevPatAtoms   = NPatAtom;
      PrevStrucAtoms = NStrucAtom;
         
#ifdef DEBUG
      fprintf(stderr, "\nPattern atoms are:\n");
      for(j=0;j<NPatAtom;j++)
      {
         fprintf(stderr, "Property: %s; Distance: %s\n",
                 PatAtom[j].property, PatAtom[j].dist);
      }

      fprintf(stderr, "\nStructure atoms are:\n");
      for(j=0;j<NStrucAtom;j++)
      {
         fprintf(stderr, "Property: %s; Distance: %s\n",
                 StrucAtom[j].property, StrucAtom[j].dist);
      }
#endif      

      /* For each atom in the structure look to see if the bit string is
         not found in the pattern. If not found, kill the atom
      */
      for(j=0; j<NStrucAtom; j++)
      {
         BOOL Found = FALSE;

         for(k=0; k<NPatAtom; k++)
         {
	    if(!strncmp(StrucAtom[j].property,
                        PatAtom[k].property, MAXPROP))
	    {
	       if(!Compare(StrucAtom[j].dist,PatAtom[k].dist,MAXDIST))
	       {
		  Found = TRUE;
		  break;
	       }
	    }
	 }

         if(!Found)
	 {
            KillAtom(StrucAtom[j].resid, struc, nstruc);
	    fprintf(stderr, "Structure atom: %s %s killed\n",
                    StrucAtom[j].resnam, StrucAtom[j].resid);
         }
      }
   }

   PrintResults(out, NPatAtom, PatAtom, NStrucAtom, StrucAtom);
}

/************************************************************************/
/*>void KillAtom(char *resid, DATA *struc, int nstruc)
   ----------------------------------------------------
   Set the dead flag in the data array for an atom;

   19.11.93 Original   By: ACRM
*/
void KillAtom(char *resid, DATA *data, int ndata)
{
   int i;

   for(i=0; i<ndata; i++)
   {
      if(!strcmp(resid, data[i].resid[0]) ||
         !strcmp(resid, data[i].resid[1]))
      {
         data[i].dead = TRUE;
      }
   }
}

/************************************************************************/
/*>void TrimBitStrings(int npat, ATOM *pat, int nstruc, ATOM *struc)
   -----------------------------------------------------------------
   Search the bit strings of the pattern and remove any distance flags
   which never occur in the structure and vice versa

   19.11.93 Original   By: ACRM
*/
void TrimBitStrings(int npat, ATOM *pat, int nstruc, ATOM *struc)
{
   int  i, j;
   BOOL PatHit,
        StrucHit;

   /* For each distance in the distance bit string                      */
   for(i=0; i<MAXDIST; i++)
   {
      PatHit = StrucHit = FALSE;

      /* Search the pattern array for this distance being flagged       */
      for(j=0; j<npat; j++)
      {
         if(pat[j].dist[i] == '1')
	 {
            PatHit = TRUE;
            break;
         }
      }

      /* Search the structure array for this distance being flagged     */
      for(j=0; j<nstruc; j++)
      {
         if(struc[j].dist[i] == '1')
	 {
            StrucHit = TRUE;
            break;
         }
      }

      /* If there was a hit in pattern, but not structure, kill refs in
         pattern
      */
      if(PatHit && !StrucHit)
      {
         for(j=0; j<npat; j++)
	 {
            pat[j].dist[i] = '0';
         }
      }

      /* If there was a hit in structure, but not pattern, kill refs in
         structure
      */
      if(StrucHit && !PatHit)
      {
         for(j=0; j<nstruc; j++)
	 {
            struc[j].dist[i] = '0';
         }
      }
   }
}

/************************************************************************/
/*>void PrintResults(FILE *out,
                     int NPatAtom,   ATOM *PatAtom, 
                     int NStrucAtom, ATOM *StrucAtom)
   --------------------------------------------------
   Run through the pattern atoms and, for each, print the best match 
   from the structure atoms

   22.11.93 Original   By: ACRM
*/
void PrintResults(FILE *out,
                  int NPatAtom,   ATOM *PatAtom, 
                  int NStrucAtom, ATOM *StrucAtom)
{
   int i;

   for(i=0; i<NPatAtom; i++)
      PrintBestMatch(out, PatAtom, i, StrucAtom, NStrucAtom);
}

/************************************************************************/
/*>void PrintBestMatch(FILE *out,
                       ATOM *PatAtom,   int PatIndex, 
                       ATOM *StrucAtom, int NStrucAtom)
   ----------------------------------------------------
   Print the best macth from the structure for this pattern atom

   22.11.93 Original   By: ACRM
*/
void PrintBestMatch(FILE *out,
                    ATOM *PatAtom,   int PatIndex, 
                    ATOM *StrucAtom, int NStrucAtom)
{
   int  j,
        best = -1;
   REAL score, 
        BestScore = 0.0;

   for(j=0; j<NStrucAtom; j++)
   {
      if(!strncmp(PatAtom[PatIndex].property,
                  StrucAtom[j].property, MAXPROP) &&
         !Compare(PatAtom[PatIndex].dist,
                  StrucAtom[j].dist,     MAXDIST))
      {
         score = CalcScore(PatAtom[PatIndex].dist, StrucAtom[j].dist,
                           MAXDIST);
         if(score > BestScore)
         {
            BestScore = score;
            best = j;
         }
      }
   }

   /* If we got a best score, print it out                              */
   if(best != (-1))
   {
      fprintf(out, "Pattern: %s %s matches Structure: %s %s\n",
              PatAtom[PatIndex].resnam, PatAtom[PatIndex].resid,
              StrucAtom[best].resnam, StrucAtom[best].resid);
   }
}

/************************************************************************/
/*>BOOL Compare(char *str1, chsr *str2, int len)
   ---------------------------------------------
   Compares two bitstrings (1/0 character strings) using a 75% identity
   requirement

   22.11.93 Original   By: ACRM
*/
BOOL Compare(char *str1, char *str2, int len)
{
   int  i,
        CountStr1 = 0,
        CountStr2 = 0,
        NMatch    = 0;
   REAL Accuracy;

   for(i=0; i<len; i++)
   {
      if(str1[i] == '1') CountStr1++;
      if(str2[i] == '1') CountStr2++;

      if((str1[i] == '1') && (str2[i] == '1')) NMatch++;
   }

   if((Accuracy = ((REAL)100.0 * (REAL)NMatch / 
                   (REAL)MAX(CountStr1, CountStr2))) >= 
      gAccuracy)
   {
      return(FALSE);
   }

   return(TRUE);
}

/************************************************************************/
/*>REAL CalcScore(char *str1, char *str2, int len)
   -----------------------------------------------
   Calculate the score for this match. Much the same as compare, but 
   returns the percentage score rather than a BOOL

   22.11.93 Original   By: ACRM
*/
REAL CalcScore(char *str1, char *str2, int len)
{
   int i,
       CountStr1 = 0,
       CountStr2 = 0,
       NMatch    = 0;

   for(i=0; i<len; i++)
   {
      if(str1[i] == '1') CountStr1++;
      if(str2[i] == '1') CountStr2++;

      if((str1[i] == '1') && (str2[i] == '1')) NMatch++;
   }

   return((REAL)100.0 * (REAL)NMatch / (REAL)MAX(CountStr1, CountStr2));
}
