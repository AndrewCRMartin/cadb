/*************************************************************************

   Program:    makecadb
   File:       makecadb.c
   
   Version:    V1.2
   Date:       18.01.02
   Function:   Create a CA distance matrix database from a PDB directory
   
   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 1998-2002
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   This is a reimplementation of the makedb method from my thesis.

   As a rough estimate, using 20 distances, the resulting database file
   is 25% larger than the source PDB files.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  06.10.98 Original
   V1.1  11.01.02 Added check that structure contains some CA atoms
   V1.2  18.01.02 Added limit on maximum number of PDB files read

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define DEF_NDIST 20

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void ProcessAllFiles(FILE *out, char *pdbdir, int ndist, int limit);
void ProcessFile(FILE *out, char *filename, int ndist);
void CalcDistances(FILE *out, char *pdbcode, PDB **pdbidx, 
                   int natoms, int ndist);
BOOL ParseCmdLine(int argc, char **argv, char *pdbdir, char *outfile, 
                  int *ndist, int *limit);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   The main routine

   06.10.98 Original   By: ACRM
   18.01.02 Added limit
*/
int main(int argc, char **argv)
{
   FILE *out = stdout;
   char outfile[MAXBUFF],
        pdbdir[MAXBUFF];
   int  ndist,
        limit = 0;
   time_t tm;
   
   if(ParseCmdLine(argc, argv, pdbdir, outfile, &ndist, &limit))
   {
      if(OpenStdFiles(NULL, outfile, NULL, &out))
      {
         fprintf(out,"!PDBDIR %s\n",pdbdir);
         fprintf(out,"!NDIST  %d\n",ndist);
         time(&tm);
         fprintf(out,"!DATE   %s\n",ctime(&tm));
         
         ProcessAllFiles(out, pdbdir, ndist, limit);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>void ProcessAllFiles(FILE *out, char *pdbdir, int ndist, int limit)
   -------------------------------------------------------------------
   Inputs:     FILE   *out     Output file pointer
               char   *pdbdir  String pointer containing name of PDB
                               directory
               int    ndist    Number of distances to calculate
               int    limit    Max number of files to process

   Steps through each file in the specified directory calling 
   ProcessFile() on each one.

   06.10.98 Original   By: ACRM
   18.01.02 Added limit
*/
void ProcessAllFiles(FILE *out, char *pdbdir, int ndist, int limit)
{
   DIR           *dp;
   struct dirent *dent;
   char          filename[MAXBUFF];
   int           count = 0;
   
   
   if((dp=opendir(pdbdir))!=NULL)
   {
      while((dent=readdir(dp))!=NULL)
      {
         if((!limit) || (count < limit))
         {
            sprintf(filename,"%s/%s",pdbdir,dent->d_name);
            ProcessFile(out, filename, ndist);
            count++;
         }
      }
      closedir(dp);
   }
}


/************************************************************************/
/*>void ProcessFile(FILE *out, char *filename, int ndist)
   ------------------------------------------------------
   Inputs:     FILE   *out         Output file pointer
               char   *filename    PDB file to be processed
               int    ndist        Number of constraints to calculate

   Reads the specified PDB file, select out the CA atoms and call
   CalcDistances() to calculate distance constraints and write results 
   to the outfile.

   06.10.98 Original   By: ACRM
   11.01.02 Added check that SelectCaPDB() found some atoms
*/
void ProcessFile(FILE *out, char *filename, int ndist)
{
   FILE *fp;
   PDB  *pdb, **pdbidx;
   int  natoms;
   char *pdbcode;

   if((fp=fopen(filename,"r"))!=NULL)
   {
      if((pdbcode = FNam2PDB(filename))!=NULL)
      {
         if((pdb=ReadPDBAtoms(fp,&natoms))!=NULL)
         {
            if((pdb=SelectCaPDB(pdb))!=NULL)
            {
               pdbidx=IndexPDB(pdb, &natoms);
               
               CalcDistances(out, pdbcode, pdbidx, natoms, ndist);
               
               free(pdbidx);
               FREELIST(pdb, PDB);
            }
         }
      }
      fclose(fp);
   }
}


/************************************************************************/
/*>void CalcDistances(FILE *out, char *pdbcode, PDB **pdbidx, 
                      int natoms, int ndist)
   ----------------------------------------------------------
   Inputs:     FILE   *out         Output file pointer
               char   *pdbcode     PDB code derived from filename
               PDB    **pdbidx     Array of PDB pointers
               int    natoms       Number of atoms in array
               int    ndist        Number of constraints to calculate

   Calculate the ndist distances between CA atoms and write them to
   the output file.

   06.10.98 Original   By: ACRM
*/
void CalcDistances(FILE *out, char *pdbcode, PDB **pdbidx, 
                   int natoms, int ndist)
{
   int  atnum, 
        currentAtom = 0,
        i = 0,
        firstAtom = 0;
   char chain = pdbidx[0]->chain[0],
        PrintChain;
   REAL dist;

   PrintChain = ((chain==' ')?'-':chain);
   
   for(atnum=0; atnum<natoms; atnum++)
   {
      if(pdbidx[atnum]->chain[0] != chain)
      {
         chain = pdbidx[atnum]->chain[0];
         firstAtom = atnum;
         PrintChain = ((chain==' ')?'-':chain);
      }
      fprintf(out,"%4s.%c.%d%c ",
              pdbcode,
              PrintChain,
              pdbidx[atnum]->resnum,
              pdbidx[atnum]->insert[0]);

      /* Do the DP (forward) distances                                */
      for(i=1; i<=ndist; i++)
      {
         currentAtom = atnum+i;
         if((currentAtom >= natoms) ||
            (pdbidx[currentAtom]->chain[0] != chain))
         {
            dist = (-1.0);
         }
         else
         {
            dist = DIST(pdbidx[atnum], pdbidx[currentAtom]);
         }
         fprintf(out, "%.2f ", dist);
      }

      /* Do the DM (backward) distances                               */
      for(i=1; i<=ndist; i++)
      {
         currentAtom = atnum-i;
         if(currentAtom < firstAtom)
         {
            dist = (-1.0);
         }
         else
         {
            dist = DIST(pdbidx[atnum], pdbidx[currentAtom]);
         }
         fprintf(out, "%.2f ", dist);
      }
      
      fprintf(out,"\n");
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *pdbdir, char *outfile, 
                     int *ndist, int *limit)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *pdbdir      PDB directory to scan
            char   *outfile     Output file (or blank string)
            int    *ndist       Number of distances
            int    *limit       Max number of PDB files to read
   Returns: BOOL                Success?

   Parse the command line
   
   06.10.98 Original    By: ACRM
   18.01.02 Added -l
*/
BOOL ParseCmdLine(int argc, char **argv, char *pdbdir, char *outfile, 
                  int *ndist, int *limit)
{
   argc--;
   argv++;

   if(argc == 0)
      return(FALSE);

   *ndist = DEF_NDIST;
   outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            argc--;
            argv++;
            sscanf(argv[0],"%d",ndist);
            break;
         case 'l':
            argc--;
            argv++;
            sscanf(argv[0],"%d",limit);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there is 1 or 2 arguments left                 */
         if((argc < 1) || (argc > 2))
            return(FALSE);
         
         strcpy(pdbdir, argv[0]);

         argc--; argv++;
         if(argc)
         {
            strcpy(outfile, argv[0]);
         }

         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   06.10.98 Original   By: ACRM
   11.01.02 V1.1
   18.01.02 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\nmakecadb V1.2 (c) 1998-2002, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: makecadb [-d ndist] [-l limit] pdbdir \
[outfile]\n");
   fprintf(stderr,"       -d Specify number of distances (Default: %d)\n",
           DEF_NDIST);
   fprintf(stderr,"       -l Limit the maximum number of PDB files read\n");

   fprintf(stderr,"\nCreates a C-alpha distance matrix database for \
use with searchdb\n");
   fprintf(stderr,"PDB files are searched from the specified pdbdir. \
Output is to\n");
   fprintf(stderr,"standard output or to the specified outfile.\n\n");
}

