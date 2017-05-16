/*********************************************************************
 * Clustal Omage - Multiple sequence alignment
 *
 * Copyright (C) 2010 University College Dublin
 *
 * Clustal-Omega is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This file is part of Clustal-Omega.
 *
 ********************************************************************/

/*
 *  RCS $Id: main.cpp 289 2013-09-17 10:09:37Z fabian $
 */

/*
 * We are using a mix of C and C++, which means that linking has to be
 * done with a C++ compiler. By using this "fake" main c++ function,
 * automake is convinced to use a C++ compiler for linking.
 *
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "mymain.h"
#include "clustal/util.h"
#include "squid/squid.h"
}


/**
 * @brief Convert an old Clustal command line parameter in the form of
 * [-/]param[=value] to new parameter if possible
 *
 * @param[out] iNewArgC_p
 * "argc" which will be incremented for each successfully converted option
 * @param[out] ppcNewArgV_p
 * "argv" to which each successfully converted options will be added
 * (caller has to free)
 * @param[in] pcOldArg
 * The old parameter and value command line option
 *
 */
void
ConvertOldCmdLineArg(int *iNewArgC_p, char ***ppcNewArgV_p, char *pcOldArg)
{
    char *pcOldParam, *pcOldValue, *pcOldArgCopy;
    char zcNotImplementedMsg[] = "WARNING: Invalid old command line option";
    
    pcOldArgCopy = CkStrdup(pcOldArg); 
    pcOldParam = strtok(pcOldArgCopy, "=");
    pcOldValue = strtok(NULL, "=");


    /* go through all options in order of appearance in clustalw2 -help
     *
     */
     
        /* data
         *
         */
    if (STR_NC_EQ("INFILE", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-i");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

    } else if (STR_NC_EQ("PROFILE1", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--profile1");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

    } else if (STR_NC_EQ("PROFILE2", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--profile2");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

        /* verbs
         *
         */

        /* missing:
         * OPTIONS
         */
        
    } else if (STR_NC_EQ("HELP", &pcOldParam[1])
               || STR_NC_EQ("CHECK", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-h");
        
    } else if (STR_NC_EQ("FULLHELP", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-h");
        
    } else if (STR_NC_EQ("ALIGN", &pcOldParam[1])) {
        char msg[] = "WARNING: The ALIGN option is default in Clustal Omega";
        fprintf(stderr, "%s\n", msg);

        /* missing:
         * TREE
         * PIM
         * BOOTSTRAP
         * CONVERT
         */
        
        /* parameters
         *
         */
    } else if (STR_NC_EQ("INTERACTIVE", &pcOldParam[1])) {
        char msg[] = "WARNING: There is no interactive command-line menu in Clustal Omega";
        fprintf(stderr, "%s\n", msg);
        /* trigger help */
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-h");
        
    } else if (STR_NC_EQ("QUICKTREE", &pcOldParam[1])) {
        char msg[] = "WARNING: The QUICKTREE (i.e. k-tuple distance) option is default in Clustal Omega";
        fprintf(stderr, "%s\n", msg);
        
    } else if (STR_NC_EQ("TYPE", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-t");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

        /* NEGATIVE */
        
    } else if (STR_NC_EQ("OUTFILE", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-o");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

        
    } else if (STR_NC_EQ("OUTPUT", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--outfmt");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

        /* missing:
         * OUTORDER
         * CASE
         * SEQNOS
         * SEQNO_RANGE
         * RANGE
         */
                   
    } else if (STR_NC_EQ("MAXSEQLEN", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--maxseqlen");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);
        
    } else if (STR_NC_EQ("QUIET", &pcOldParam[1])) {
        char msg[] = "WARNING: The QUIET option is default in Clustal Omega";
        fprintf(stderr, "%s\n", msg);

        /* missing:
         * STATS
         */
        
        /* fast pariwise alignment
         *
         */

        /* missing:
         * KTUPLE
         * TOPDIAGS
         * WINDOW
         * PAIRGAP
         * SCORE
         */
        
        /* slow pairwise alignments 
         *
         */

        /* missing:
         * PWMATRIX
         * PWDNAMATRIX
         * PWGAPOPEN
         * PWGAPEXT
         */
        
        /* multiple alignments
         *
         */
    } else if (STR_NC_EQ("NEWTREE", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--guidetree-out");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);
        
    } else if (STR_NC_EQ("USETREE", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--guidetree-in");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);

        /* missing:
         * MATRIX
         * DNAMATRIX
         * GAPOPEN
         * GAPEXT
         * ENDGAPS
         * GAPDIST
         * NOGAP
         * NOHGAP
         * HGAPRESIDUES
         * MAXDIV
         * TYPE already handled above
         * TRANSWEIGHT
         * ITERATION
         * NUMITER
         * NOWEIGHTS
         */

        /* profile alignments
         *
         */

        /* missing:
         * PROFILE
         * NEWTREE1
         * NEWTREE2
         * USETREE1
         * USETREE2
         */
        
        /* sequence to profile alignments 
         *
         */
    } else if (STR_NC_EQ("SEQUENCES", &pcOldParam[1])) {
        fprintf(stderr, "WARNING: %s: %s\n", zcNotImplementedMsg, pcOldArg);

        /* SEQUENCES and NEWTREE already handled above */
        
        /* structure alignments
         *
         */

        /* missing:
         * NOSECSTR1
         * NOSECSTR2
         * SECSTROUT
         * HELIXGAP
         * STRANDGAP
         * LOOPGAP
         * TERMINALGAP
         * HELIXENDIN
         * HELIXENDOUT
         * STRANDENDIN
         * STRANDENDOUT
         */
        
        /* trees
         *
         */

        /* missing:
         * OUTPUTTREE
         * SEED
         * KIMURA
         * TOSSGAPS
         * BOOTLABELS
         */
#if 0        
    } else if (STR_NC_EQ("CLUSTERING", &pcOldParam[1])) {
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-c");
        if (NULL != pcOldValue)
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(pcOldValue);
#endif
        
    } else {
        fprintf(stderr,
                "WARNING: Unsupported old command line option '%s' will be ignored (may change default output stream and format)\n",
                pcOldArg);
    }
    
	/* FIXME: if not outfile was given, than the old default was to create a
	 * filename based on the input filename but with its extension replaced by
	 * aln. If the input already had the extension aln then the input was
	 * overwritten. What to do here? Strictly mimic the old behaviour?
	 */
	
    free(pcOldArgCopy);
}
/* end ConvertOldCmdLineArg */



/**
 * @brief Convert old command line usage to new one. Mix of old and
 * new style will be tolerated
 *
 * @param[out] iNewArgC_p
 * The updated "argc"
 * @param[out] ppcNewArgV_p
 * The updated "argv". Caller has to free.
 * @param[in] argc
 * Original "argc"
 * @param[in] argv
 * Original argv
 * 
 * @note old style parameters look like this:
 *  [/-]param[=value]
 * new style parameters:
 *  -p [value]
 *  --param [value]
 *
 */
void
ConvertOldCmdline(int *iNewArgC_p, char ***ppcNewArgV_p, int argc, char **argv)
{
    bool bOldCmdlineDetected = false;
    int i; /* aux */

    /* we can have at most 2*argc converted arguments, plus the few
     * that we set by default (.e.g --force)
     */
    (*ppcNewArgV_p) = (char **) CKCALLOC(argc*2+10, sizeof(char*));
    
    /* copy first arg which is program name */
    (*ppcNewArgV_p)[0] = CkStrdup(argv[0]);
    *iNewArgC_p = 1;
        
    for (i=1; i<argc; i++) {
        bool bNewStyle = false;
        
        if (strlen(argv[i])<=2) {
            /* e.g. -i (param) or just numbers (value) */
            bNewStyle = true;
            
        } else if (strlen(argv[i])>2) {
            if (argv[i][0] == '-' && argv[i][1] == '-') {
                /* new style long opts */
                bNewStyle = true;
                
            } else if (argv[i][0]=='/' && (NULL!=strchr(&argv[i][1], '/'))) {
                /* Slash used to be a valid replacement for dash in
                 * Clustal<=2, but could also be file in new style. If
                 * we find at least two slashes, one at the beginning,
                 * it should be a filename and therefore new style */
                bNewStyle = true;
                
            } else if (argv[i][0] != '/' && argv[i][0] != '-') {
                /* old style opts always start with slash or dash */
                bNewStyle = true;
                
            }
        }

        /* copy and continue if new style arg or attempt to convert
         * old style arg
         */
        if (bNewStyle) {
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup(argv[i]);
        } else {
            ConvertOldCmdLineArg(iNewArgC_p, ppcNewArgV_p, argv[i]);
            /*LOG_DEBUG("old command line arg: %s", argv[i]);*/
            bOldCmdlineDetected = true;
        }
    }
    
    if (bOldCmdlineDetected) {
        bool bOutfileOptSet = FALSE;
        bool bOutFormatOptSet = FALSE;
        
        
        /* old clustal used to write to a file called in.aln by
         *  default. set if not
         * explicitely requested otherwisee
         */
        for (i=0; i<*iNewArgC_p; i++) {
            const char *pcOpt = "-o";
            if (strlen(pcOpt) <= strlen((*ppcNewArgV_p)[i])) {
                if (0 == strncmp((*ppcNewArgV_p)[i], pcOpt, strlen(pcOpt))) {
                    bOutfileOptSet = TRUE;
                    break;
                }
            }
        }
        if (FALSE == bOutfileOptSet) {
#ifdef TOO_LAZY_TO_IMPLEMENT_JUST_USING_DEFAULT_NAME_INSTEAD
            char *pcDotPos = NULL;
            char *pcInfileOpt = NULL;

            /* get infile arg and find last dot in it. if found replace
             * everything after with "aln", otherwise just add "aln"
             */
            for (i=0; i<*iNewArgC_p; i++) {
                const char *pcOpt = "-i";
                if (strlen(pcOpt) <= strlen((*ppcNewArgV_p)[i])) {
                    if (0 == strncmp((*ppcNewArgV_p)[i], pcOpt, strlen(pcOpt))) {
                        if (*iNewArgC_p<= i+1) {
                            fprintf(stderr,
                                    "Oups...error while trying to convert old commandline (%s).\n",
                                    "No more arguments left after -i");
                            exit(1);
                        }
                        pcInfileOpt = (*ppcNewArgV_p)[i+1];
                        break;
                    }
                }
            }
            if (NULL == pcInfileOpt) {
                fprintf(stderr,
                        "Oups...error while trying to convert old commandline (%s)\n",
                        "No infile opt found");
                exit(1);
            }

            fprintf(stderr, "FIXME: unfinished\n");
            exit(1);
#endif
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-o");
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("clustal.aln");
        }

        
        /* old clustal used the clustal format by default. set if not
         * explicitely requested otherwisee
         */
        for (i=0; i<*iNewArgC_p; i++) {
            const char *pcOpt = "--outfmt";
            if (strlen(pcOpt) <= strlen((*ppcNewArgV_p)[i])) {
                if (0 == strncmp((*ppcNewArgV_p)[i], pcOpt, strlen(pcOpt))) {
                    bOutFormatOptSet = TRUE;
                    break;
                }
            }
        }
        if (FALSE == bOutFormatOptSet) {
            (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--outfmt=clustal");
        }

        
        /* old clustal was verbose by default
         */
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("-v");
        
        /* old clustal used to overwrite files by default
         */
        (*ppcNewArgV_p)[(*iNewArgC_p)++] = CkStrdup("--force");
        
        fprintf(stderr,
                "WARNING: Your old-style command-line options were converted to: ");
        for (i=0; i<*iNewArgC_p; i++) {
            fprintf(stderr, " %s", (*ppcNewArgV_p)[i]);
        }
        fprintf(stderr, "\n");
    }

}
/* end ConvertOldCmdline */



int
main(int argc, char **argv)
{
    int new_argc = 0;
    char **new_argv = NULL;
    int i; /* aux */
    
    ConvertOldCmdline(&new_argc, &new_argv, argc, argv);
    
    MyMain(new_argc, new_argv);
    
    for (i=0; i<new_argc; i++) {
        free(new_argv[i]);
    }
    free(new_argv);
    
    return 0;
}
