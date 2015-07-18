/*************************************************************************

   Program:    searchcadb
   File:       searchcadb.c
   
   Version:    V1.0
   Date:       08.10.98
   Function:   Search a CA distance matrix database
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1998
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
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
   This is a reimplementation of the searchdb method from my thesis.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0 08.10.98 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gdbm.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/parse.h"
#include "bioplib/general.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/* Defines for Keyword parser                                           */
#define KEY_DATABASE 0
#define KEY_DP       1
#define KEY_DM       2
#define KEY_END      3
#define KEY_LENGTH   4
#define KEY_QUIT     5
#define KEY_HELP     6
#define NCOMM        7
#define MAXSTRPARAM  1
#define MAXREALPARAM 3

/* Structure to store distance constraints                              */
typedef struct _constraint
{
   struct _constraint *next;
   REAL min, max;
   int cons;
}  CONSTRAINT;

/* Make a key to store in a dbm hash                                    */
#define MAKEDBMKEY(MK_datum, MK_key) \
        do { (MK_datum).dptr = (MK_key); \
             (MK_datum).dsize = 1+strlen(MK_key); } while(0)
        

/************************ The ERRPROMPT macro ***************************/
/* Default is just to print a string as a prompt                        */
#define ERRPROMPT(in,x) fprintf(stderr,"%s",(x))

/* More intelligent prompts for systems where we know the FILE structure*/
#ifdef __sgi
#  undef ERRPROMPT
#  define ERRPROMPT(in,x) do{if(isatty((in)->_file)) \
                         fprintf(stderr,"%s",(x));}while(0)
#endif
#ifdef __linux__
#  undef ERRPROMPT
#  define ERRPROMPT(in,x) do{if(isatty((in)->_fileno)) \
                         fprintf(stderr,"%s",(x));}while(0)
#endif


/************************************************************************/
/* Globals
*/
KeyWd      gKeys[NCOMM];
char       *gStrParam[MAXSTRPARAM];
REAL       gRealParam[MAXREALPARAM];
CONSTRAINT *gPosConsList = NULL,
           *gNegConsList = NULL;
DBM        *gDbm = NULL;
int        gLoopLength = 0;


/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *InFile, char *OutFile);
BOOL SetupParser(void);
BOOL ParseInputFile(FILE *in, FILE *out);
BOOL StorePosConstraint(int cons, REAL mindist, REAL maxdist);
BOOL StoreNegConstraint(int cons, REAL mindist, REAL maxdist);
BOOL RunSearch(FILE *DBfp, int ndist, FILE *out);
void ReadArrayFromBuffer(char *buffer, int ndist, REAL *distArray);
BOOL RecordOK(REAL *distArray, int ndist, CONSTRAINT *ConsList);
BOOL InSameChain(char *currentKey, char *prevKey);
void FlagPosOK(char *currentKey);
void FlagNegBad(char *prevKey);
void DisplayResults(FILE *out);
void ShowHelp(void);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   08.10.98 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF];
   BOOL Success;
   
   if(ParseCmdLine(argc, argv, InFile, OutFile))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if(SetupParser())
         {
            /* initialise a DBM file                                   */
            if((gDbm = dbm_open("/tmp/scadb.dbm", 
                                O_RDWR|O_CREAT|O_TRUNC, 
                                S_IRWXU)) == NULL)
            {
               fprintf(stderr,"Can't open NDBM file for writing\n");
               return(1);
            }
            

            Success = ParseInputFile(in,out);

            dbm_close(gDbm);
            unlink("/tmp/scadb.dbm");
            
            if(!Success)
               return(1);
         }
         else
         {
            fprintf(stderr,"No memory for parser strings\n");
            return(1);
         }
      }
   }
   else
   {
      Usage();
   }

   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *InFile, char *OutFile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *InFile      Input file (or blank string)
            char   *OutFile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   08.10.98 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *InFile, char *OutFile)
{
   argc--;
   argv++;

   InFile[0] = '\0';
   OutFile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1 or 2 arguments left                  */
         if((argc < 1) || (argc > 2))
            return(FALSE);
         
         strcpy(InFile, argv[0]);

         argc--; argv++;
         if(argc)
         {
            strcpy(OutFile, argv[0]);
         }

         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL SetupParser(void)
   ----------------------
   Returns:    BOOL           Success?

   Sets up the command parser.

   08.10.98 Original   By: ACRM
*/
BOOL SetupParser(void)
{
   int i;
   
   for(i=0; i<MAXSTRPARAM; i++)
   {
      if((gStrParam[i] = (char *)malloc(MAXBUFF * sizeof(char)))==NULL)
      {
         return(FALSE);
      }
   }
   
   MAKEKEY(gKeys[KEY_DATABASE], "DATABASE", STRING, 1);
   MAKEKEY(gKeys[KEY_DP],       "DP",       NUMBER, 3);
   MAKEKEY(gKeys[KEY_DM],       "DM",       NUMBER, 3);
   MAKEKEY(gKeys[KEY_END],      "END",      NUMBER, 0);
   MAKEKEY(gKeys[KEY_LENGTH],   "LENGTH",   NUMBER, 1);
   MAKEKEY(gKeys[KEY_QUIT],     "QUIT",     NUMBER, 0);
   MAKEKEY(gKeys[KEY_HELP],     "HELP",     NUMBER, 0);

   return(TRUE);
}


/************************************************************************/
/*>BOOL ParseInputFile(FILE *in, FILE *out)
   ----------------------------------------
   Inputs:     FILE   *in       Input control file
   Outputs:    FILE   *out      Results file
   Returns:    BOOL             Success?

   Runs through the control file, handling specified commands and 
   calling routines to act on them.

   08.10.98 Original   By: ACRM
*/
BOOL ParseInputFile(FILE *in, FILE *out)
{
   char buffer[MAXBUFF];
   FILE *DBfp = NULL;
   int  ndist = 20,
        i;
   
   ERRPROMPT(in,"SEARCHCADB> ");
   
   while(fgets(buffer,MAXBUFF,in))
   {
      TERMINATE(buffer);

      switch(parse(buffer,NCOMM,gKeys,gRealParam,gStrParam))
      {
      case PARSE_ERRC:
         fprintf(stderr,"Error in command: %s\n",buffer);
         break;
      case PARSE_ERRP:
         fprintf(stderr,"Error in parameters: %s\n",buffer);
         break;
      case KEY_DATABASE:
         if(DBfp != NULL)
         {
            fprintf(stderr,"Database already open, command ignored\n");
         }
         else
         {
            if((DBfp=fopen(gStrParam[0],"r"))==NULL)
            {
               fprintf(stderr,"Can't open database: %s\n",gStrParam[0]);
            }
            else
            {
               for(i=0; i<3; i++)
               {
                  if(fgets(buffer,MAXBUFF,DBfp))
                  {
                     TERMINATE(buffer);
                     if(!strncmp(buffer,"!NDIST",6))
                     {
                        sscanf(buffer+6,"%d",&ndist);
                        break;
                     }
                  }
               }
            }
         }
         break;
      case KEY_DP:
         if(!StorePosConstraint((int)gRealParam[0],
                                gRealParam[1],gRealParam[2]))
         {
            fprintf(stderr,"No memory for constraint list\n");
            return(FALSE);
         }
         break;
      case KEY_DM:
         if(!StoreNegConstraint((int)gRealParam[0],
                                gRealParam[1],gRealParam[2]))
         {
            fprintf(stderr,"No memory for constraint list\n");
            return(FALSE);
         }
         break;
      case KEY_END:
         if(gLoopLength == 0)
         {
            fprintf(stderr,"You must specify a loop length first!\n");
         }
         else
         {
            if(DBfp!=NULL)
            {
               return(RunSearch(DBfp,ndist,out));
            }
            else
            {
               fprintf(stderr,"Database must be opened first!\n");
            }
         }
         break;
      case KEY_LENGTH:
         gLoopLength = (int)gRealParam[0];
         break;
      case KEY_QUIT:
         return(TRUE);
         break;
      case KEY_HELP:
         ShowHelp();
         break;
      default:
         break;
      }
      ERRPROMPT(in,"SEARCHCADB> ");
   }
   return(TRUE);
}


/************************************************************************/
/*>BOOL StorePosConstraint(int cons, REAL mindist, REAL maxdist)
   -------------------------------------------------------------
   Inputs:     int        cons          Constraint number
               REAL       mindist       Minimum distance
               REAL       maxdist       Maximum distance
   Returns:    BOOL                     Success?
   Globals:    CONSTRAINT gPosConsList  +ve constraints linked list

   Store a positive distance constraint in the linked list.

   08.10.98 Original   By: ACRM
*/
BOOL StorePosConstraint(int cons, REAL mindist, REAL maxdist)
{
   static CONSTRAINT *c;
   
   if(gPosConsList==NULL)
   {
      INIT(gPosConsList, CONSTRAINT);
      c = gPosConsList;
   }
   else
   {
      ALLOCNEXT(c, CONSTRAINT);
   }

   if(c==NULL)
   {
      if(gPosConsList)
         FREELIST(gPosConsList, CONSTRAINT);
      return(FALSE);
   }
   
   c->cons = cons;
   c->min  = mindist;
   c->max  = maxdist;
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL StoreNegConstraint(int cons, REAL mindist, REAL maxdist)
   -------------------------------------------------------------
   Inputs:     int        cons          Constraint number
               REAL       mindist       Minimum distance
               REAL       maxdist       Maximum distance
   Returns:    BOOL                     Success?
   Globals:    CONSTRAINT gNegConsList  -ve constraints linked list

   Store a negative distance constraint in the linked list.

   08.10.98 Original   By: ACRM
*/
BOOL StoreNegConstraint(int cons, REAL mindist, REAL maxdist)
{
   static CONSTRAINT *c;
   
   if(gNegConsList==NULL)
   {
      INIT(gNegConsList, CONSTRAINT);
      c = gNegConsList;
   }
   else
   {
      ALLOCNEXT(c, CONSTRAINT);
   }

   if(c==NULL)
   {
      if(gNegConsList)
         FREELIST(gNegConsList, CONSTRAINT);
      return(FALSE);
   }
   
   c->cons = cons;
   c->min  = mindist;
   c->max  = maxdist;
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL RunSearch(FILE *DBfp, int ndist, FILE *out)
   ------------------------------------------------
   Inputs:     FILE    *DBfp       Database file pointer
               int     ndist       Number of distance constraints in file
               FILE    *out        Output file pointer
   Returns:    BOOL                Success?

   Actually runs the search. 

   Checking DP constraints is easy!

   For DM constraints we need to update the N-LoopLength record. This we
   do by keeping a cyclic list (prevKeys) of the LoopLength previous keys.
   keyPos points to the next position in which we will insert a key; once
   the list has cycled once, it is also the position of the N-LoopLength
   key. 

   We depend on the fact that the main database file contains records in 
   the correct order of the atoms!

   Allocates memory for the prevKeys array and distances array. Steps 
   through the database file parsing in the set of distances. Checks
   if the +ve constraints are OK and if so stashes this identifier in
   a DBM hash. Once we've got enough records, then checks the -ve
   constraints and if these fail, finds the beginning of the loop from
   the prevKeys cyclic array and then deletes this key from the DBM hash.

   08.10.98 Original   By: ACRM
*/
BOOL RunSearch(FILE *DBfp, int ndist, FILE *out)
{
   char *buffer    = NULL,
        **prevKeys = NULL,
        currentKey[16];
   int  bufferSize,
        keyPos     = 0;
   REAL *distArray = NULL;
   BOOL Cycled     = FALSE;

   /* Allocate string array to store the cycle of previous keys         */
   if((prevKeys = (char**)Array2D(sizeof(char),gLoopLength,16))==NULL)
   {
      fprintf(stderr,"No memory for previous key array\n");
      return(FALSE);
   }
   
   /* buffer is used to store lines read from the database file         */
   bufferSize = (2 * ndist * 7) + 100;
   if((buffer=(char *)malloc(bufferSize * sizeof(char)))==NULL)
   {
      FreeArray2D(prevKeys,gLoopLength,16);
      fprintf(stderr,"No memory for buffer to read database file\n");
      return(FALSE);
   }

   /* distArray stores the distances parsed oyt of the database file    */
   if((distArray=(REAL *)malloc(2*ndist*sizeof(REAL)))==NULL)
   {
      FreeArray2D(prevKeys,gLoopLength,16);
      fprintf(stderr,"No memory for distance array\n");
      free(buffer);
      return(FALSE);
   }

   while(fgets(buffer,bufferSize,DBfp))
   {
      TERMINATE(buffer);
      if((buffer[0] == '!') ||
         (buffer[0] == '#') ||
         !strlen(buffer))
         continue;

      sscanf(buffer,"%s",currentKey);
      strcpy(prevKeys[keyPos],currentKey);
      if(++keyPos >= gLoopLength)
      {
         keyPos = 0;
         Cycled = TRUE;
      }

      ReadArrayFromBuffer(buffer,ndist,distArray);

      if(RecordOK(distArray, 0, gPosConsList))
      {
         FlagPosOK(currentKey);
      }
      if(Cycled)
      {
         if(InSameChain(currentKey, prevKeys[keyPos]))
         {
            if(!RecordOK(distArray, ndist, gNegConsList))
            {
               FlagNegBad(prevKeys[keyPos]);
            }
         }
      }
   }

   /* Display the flagged records                                       */
   DisplayResults(out);

   free(buffer);
   return(TRUE);
}


/************************************************************************/
/*>void ReadArrayFromBuffer(char *buffer, int ndist, REAL *distArray)
   ------------------------------------------------------------------
   Inputs:     char  *buffer      Buffer read from database file
               int   ndist        Number of distances in file
   Outputs:    REAL  *distArray   Array of parsed distances

   Parses a set of distances out of the buffer into the distArray

   08.10.98 Original   By: ACRM
*/
void ReadArrayFromBuffer(char *buffer, int ndist, REAL *distArray)
{
   char *chp,
        word[16];
   int  i;
   
   chp = buffer;
   /* Junk the first word which is the identifier              */
   chp=GetWord(chp,word);
   
   /* Put the others into the distance array                   */
   i=0;
   while(((chp=GetWord(chp,word))!=NULL) && (i<ndist*2))
   {
      sscanf(word,"%lf",&(distArray[i++]));
   }
}


/************************************************************************/
/*>BOOL RecordOK(REAL *distArray, int ndist, CONSTRAINT *ConsList)
   ---------------------------------------------------------------
   Inputs:     REAL       *distArray   Set of distances for this record
               int        ndist        Number of distances
               CONSTRAINT *ConsList    Linked list of constraints
   Returns:    BOOL                    Matches constraints?

   Does a test to see if the constraints in the linked list are all
   satisfied. ndist is used as an offset when testing the -ve
   distances. It is the number of +ve distances and therefore the
   number of columns which must be skipped. i.e. set it to 0 for +ve
   constraints and to the real ndist for the -ve constraints.

   08.10.98 Original   By: ACRM
*/
BOOL RecordOK(REAL *distArray, int ndist, CONSTRAINT *ConsList)
{
   CONSTRAINT *c;

   for(c=ConsList; c!=NULL; NEXT(c))
   {
      if((distArray[c->cons + ndist - 1] < c->min) ||
         (distArray[c->cons + ndist - 1] > c->max))
      {
         return(FALSE);
      }
   }
   return(TRUE);
}

/************************************************************************/
/*>BOOL InSameChain(char *currentKey, char *prevKey)
   -------------------------------------------------
   Inputs:     char  *currentKey    Current identifier
   Outputs:    char  *prevKey       Previous identifier
   Returns:    BOOL                 In same chain?

   Tests whether 2 identifiers are in the same protein chain

   08.10.98 Original   By: ACRM
*/
BOOL InSameChain(char *currentKey, char *prevKey)
{
   if(!strncmp(currentKey, prevKey, 6))
      return(TRUE);

   return(FALSE);
}

/************************************************************************/
/*>void FlagPosOK(char *currentKey)
   --------------------------------
   Inputs:     char   *currentKey    The current resdiue identifier

   Stores the key in the DBM hash.

   08.10.98 Original   By: ACRM
*/
void FlagPosOK(char *currentKey)
{
   datum dtm,
         dtm_ok;
   MAKEDBMKEY(dtm,    currentKey);
   MAKEDBMKEY(dtm_ok, "OK");

   D("Stored key: ");
   D(currentKey);
   D("\n");
   
   dbm_store(gDbm, dtm, dtm_ok, DBM_INSERT);
}

/************************************************************************/
/*>void FlagNegBad(char *prevKey)
   ------------------------------
   Inputs:     char   *prevKey       The previous residue identifier

   Removes the key from the DBM hash

   08.10.98 Original   By: ACRM
*/
void FlagNegBad(char *prevKey)
{
   datum dtm;
   MAKEDBMKEY(dtm, prevKey);

   D("Removed key: ");
   D(prevKey);
   D("\n");

   dbm_delete(gDbm, dtm);
}

/************************************************************************/
/*>void DisplayResults(FILE *out)
   ------------------------------
   Inputs:     FILE    *out         Output file to write to

   Display the final results.
   Simply steps through the DBM hash and prints out the keys

   08.10.98 Original   By: ACRM
*/
void DisplayResults(FILE *out)
{
   datum dtm;
   
   dtm=dbm_firstkey(gDbm);
   if(dtm.dptr != NULL)
      fprintf(out,"%s\n",dtm.dptr);

   for(;;)
   {
      dtm=dbm_nextkey(gDbm);
      if(dtm.dptr == NULL)
         break;

      fprintf(out,"%s\n",dtm.dptr);
   }
}

/************************************************************************/
/*>void ShowHelp(void)
   -------------------
   Print a help message when running the program.

   08.10.98 Original   By: ACRM
*/
void ShowHelp(void)
{
   fprintf(stderr,"DATABASE dbname     Specify the database written by \
makecadb\n");
   fprintf(stderr,"LENGTH length       Specify loop length\n");
   fprintf(stderr,"DP n min max        Distance constraint from Nter of \
loop\n");
   fprintf(stderr,"DM n min max        Distance constraint from Cter of \
loop\n");
   fprintf(stderr,"END                 Run the search\n");
   fprintf(stderr,"QUIT                Exit without running the \
search\n");
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Print a usage message

   08.10.98 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nsearchcadb V1.0 (c) 1998, UCL, Dr. Andrew C.R. \
Martin\n");

   fprintf(stderr,"\nUsage: searchdb [infile [outfile]]\n");

   fprintf(stderr,"\nPerforms a search for loop conformations using \
the method of \n");
   fprintf(stderr,"Martin et al. PNAS 86(1989),9269-9272.\n");

   fprintf(stderr,"\nUsage of the program is keyword driven. Run \
searchdb and then issue \n");
fprintf(stderr,"the 'help' command for information on the available \
keywords.\n\n");
}

