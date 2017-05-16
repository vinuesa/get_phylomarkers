/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

/*********************************************************************
 * Clustal Omega - Multiple sequence alignment
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
 *  RCS $Id: mymain.c 296 2014-10-07 12:15:41Z fabian $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <argtable2.h>
#include <ctype.h>
#include <limits.h>
#include <libgen.h> /* for basename only */

/* clustal */
#include "clustal-omega.h"

#include "mymain.h"

typedef struct {

    /* Sequence input
     */
    /** sequence type (from cmdline arg) */
    int iSeqType;
    /** sequence input file. not directly used by Align() */
    char *pcSeqInfile;
    /** Dealign sequences on input. Otherwise we use the alignment
     * info for background-HMM creation */
    bool bDealignInputSeqs;

    /** Sequence input format
     */
    int iSeqInFormat;

    /* profiles: : pre-aligned sequences, whose alignment will not be changed 
     */
    /** profile 1: pre-aligned sequence input. not directly used by Align() */
    char *pcProfile1Infile ;
    /** profile 2: pre-aligned sequence input. not directly used by Align() */
    char *pcProfile2Infile;        
    /** profiles that contain no gaps are rejected, force them */
    bool bIsProfile; 
    /** up to version 1.1.1 Kimura distance was default, change default, make Kimura optional */
    /*bool bUseKimura; */
    /** distance matrix output is default, allow %-ID output */
    bool bPercentID; 

    /** Input limitations
     */
    /** maximum allowed number of input sequences */
    int iMaxNumSeq;
    /** maximum allowed input sequence length */
    int iMaxSeqLen;

    /* Alignment output
     */
    /** alignment output file */
    char *pcAlnOutfile;
    /** alignment output format */
    int iAlnOutFormat;
    /** force overwriting of files */
    bool bForceFileOverwrite;

    /* line wrapping, FS, r274 -> */
    int iWrap;
    /* residue number in Clustal format, FS, r274 -> */
    bool bResno;
    /* output order, FS, r274 -> */
    int iOutputOrder; 

    /* multithreading
     */
    /** number of threads */
    int iThreads;

    /* pseudo-count file
     */
    char *pcPseudoFile; 
    /* logging 
     */
    char *pcLogFile;

    opts_t aln_opts;

    /* changes here will have to be reflected in SetDefaultUserOpts(),
     * FreeUserOpts(), PrintUserOpts() and UserOptsLogicCheck() etc
     */
} cmdline_opts_t;



/* log-file used for non-essential logging in prLog */
FILE *prLogFile = NULL;

const char *CITATION = " Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG."
    "\n"
    " Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega."
    "\n"
    " Mol Syst Biol. 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75. PMID: 21988835.";



/**
 * @brief Sets default user/commandline options
 *
 * @param[out] opts
 * The option structure to initialise
 *
 */
void
SetDefaultUserOpts(cmdline_opts_t *opts)
{
    assert(NULL != opts);

    opts->iSeqType = SEQTYPE_UNKNOWN;
    opts->pcSeqInfile = NULL;
    opts->bDealignInputSeqs = FALSE;        
    opts->pcProfile1Infile = NULL;
    opts->pcProfile2Infile = NULL;
    opts->bIsProfile = FALSE;
    opts->aln_opts.bUseKimura = FALSE;
	opts->aln_opts.bPercID = FALSE;
    opts->aln_opts.pcHMMBatch = NULL;

    opts->iMaxNumSeq = INT_MAX;
    opts->iMaxSeqLen = INT_MAX;
    
    opts->pcAlnOutfile = NULL;
    opts->iAlnOutFormat =  MSAFILE_A2M;
    opts->bForceFileOverwrite = FALSE;
    opts->iWrap = 60;
    opts->bResno = FALSE;
    opts->iOutputOrder = INPUT_ORDER;

#ifdef HAVE_OPENMP
    /* defaults to # of CPUs */
    opts->iThreads = omp_get_max_threads();
#else
    opts->iThreads = 1;
#endif
    
    opts->pcPseudoFile = NULL;
    opts->pcLogFile = NULL;

    SetDefaultAlnOpts(& opts->aln_opts);
}
/* end of SetDefaultUserOpts() */




/**
 * @brief FIXME add doc
 *
 */
void
PrintUserOpts(FILE *prFile, cmdline_opts_t *opts) {
    
    /* keep in same order as in struct. FIXME could this be derived from argtable?
     */
    fprintf(prFile, "seq-type = %s\n", SeqTypeToStr(opts->iSeqType));
    fprintf(prFile, "seq-in-fmt = %s\n", SeqfileFormat2String(opts->iSeqInFormat));
    fprintf(prFile, "option: seq-in = %s\n", 
            NULL != opts->pcSeqInfile? opts->pcSeqInfile: "(null)");
    fprintf(prFile, "option: dealign = %d\n", opts->bDealignInputSeqs);
    fprintf(prFile, "option: profile1 = %s\n", 
            NULL != opts->pcProfile1Infile? opts->pcProfile1Infile: "(null)");
    fprintf(prFile, "option: profile2 = %s\n",
            NULL != opts->pcProfile2Infile? opts->pcProfile2Infile: "(null)");
    /*fprintf(prFile, "option: percent-id = %d\n", opts->aln_opts.bPercID);*/
    fprintf(prFile, "option: is-profile = %d\n", opts->bIsProfile);
    /*fprintf(prFile, "option: use-kimura = %d\n", opts->aln_opts.bUseKimura);*/

    fprintf(prFile, "option: max-num-seq = %d\n", opts->iMaxNumSeq);
    fprintf(prFile, "option: max-seq-len = %d\n", opts->iMaxSeqLen);
    fprintf(prFile, "option: aln-out-file = %s\n", 
            NULL != opts->pcAlnOutfile? opts->pcAlnOutfile: "(null)");
    fprintf(prFile, "option: aln-out-format = %s\n", SeqfileFormat2String(opts->iAlnOutFormat));
    fprintf(prFile, "option: force-file-overwrite = %d\n", opts->bForceFileOverwrite);
    fprintf(prFile, "option: line wrap = %d\n", opts->iWrap);
    fprintf(prFile, "option: print residue numbers = %d\n", opts->bResno);
    fprintf(prFile, "option: order alignment like input/tree = %d\n", opts->iOutputOrder);

    fprintf(prFile, "option: threads = %d\n", opts->iThreads);
    fprintf(prFile, "option: PseudoFile = %s\n", opts->pcPseudoFile);
    fprintf(prFile, "option: logFile = %s\n", opts->pcLogFile);
}
/* end of PrintUserOpts */



/**
 * @brief Frees user opt members allocated during parsing
 *
 * @param[out] user_opts
 * user options whose members are to free
 *
 * @see ParseCommandLine()
 *
 */    
void
FreeUserOpts(cmdline_opts_t *user_opts)
{

    if (NULL != user_opts->pcSeqInfile) {
        CKFREE(user_opts->pcSeqInfile);
    }
    if (NULL != user_opts->pcProfile1Infile) {
        CKFREE(user_opts->pcProfile1Infile);
    }
    if (NULL != user_opts->pcProfile2Infile) {
        CKFREE(user_opts->pcProfile2Infile);
    }
    if (NULL != user_opts->pcAlnOutfile) {
        CKFREE(user_opts->pcAlnOutfile);
    }
    if (NULL != user_opts->pcLogFile) {
        CKFREE(user_opts->pcLogFile);
    }
    if (NULL != user_opts->pcPseudoFile) {
        CKFREE(user_opts->pcPseudoFile);
    }

    FreeAlnOpts(& user_opts->aln_opts);

    return; 
}
/* end of FreeUserOpts() */




/**
 * @brief Do quick&dirty logic check of used options and call Log(&rLog, LOG_FATAL, ) in case
 * of any inconsistencies
 *
 * @param[in] opts
 * option structure to check
 *
 */
void
UserOptsLogicCheck(cmdline_opts_t *opts)
{
    /* sequence input
     *
     */
    if (NULL == opts->pcSeqInfile) {
        if (NULL == opts->pcProfile1Infile && NULL == opts->pcProfile2Infile) {
            Log(&rLog, LOG_FATAL, "No sequence input was provided. For more information try: --help");
        }
    } else {
        if (NULL != opts->pcProfile1Infile && NULL != opts->pcProfile2Infile) {
            Log(&rLog, LOG_FATAL, "Can't align two profile alignments AND a 'normal' sequence file");
        }
    }
    /* if a profile was given it should always be no 1, not 2 */
    if (NULL == opts->pcProfile1Infile && NULL != opts->pcProfile2Infile) {
        Log(&rLog, LOG_FATAL, "Got a second profile, but no first one.");
    }

    /* alignment output
     */
    if (rLog.iLogLevelEnabled < LOG_WARN && NULL==opts->pcAlnOutfile && NULL==opts->pcLogFile) {
        Log(&rLog, LOG_FATAL, "%s %s",
              "You requested alignment output to stdout and verbose logging.",
              " Alignment and log messages would get mixed up.");
    }
    /* distance matrix output impossible in mBed mode, only have clusters, FS, r254 ->  */
#if 1
    /* allow distance matrix output if initial mBed but subsequent full matrix during iteration, FS, r274 -> */
    if (NULL != opts->aln_opts.pcDistmatOutfile){
        if ( (TRUE == opts->aln_opts.bUseMbed) && (opts->aln_opts.iNumIterations < 1) ){
            Log(&rLog, LOG_FATAL, "Distance Matrix output not possible in mBed mode.");
        }
        if ( (TRUE == opts->aln_opts.bUseMbed) && (TRUE == opts->aln_opts.bUseMbedForIteration) ){
            Log(&rLog, LOG_FATAL, "Distance Matrix output not possible in mBed mode.");
        }
        if ( (TRUE == opts->aln_opts.bUseMbed) && (opts->aln_opts.iNumIterations > 0) && 
             (opts->aln_opts.iMaxGuidetreeIterations < 1) ){
            Log(&rLog, LOG_FATAL, "Distance Matrix output not possible in mBed mode.");
        }
    }
#else
    if ( (NULL != opts->aln_opts.pcDistmatOutfile) && (TRUE == opts->aln_opts.bUseMbed) ) {
        Log(&rLog, LOG_FATAL, "Distance Matrix output not possible in mBed mode.");
    }
#endif

    /* percentage identity cannot be printed in Kimura mode */
    if ( (TRUE == opts->aln_opts.bUseKimura) && (TRUE == opts->aln_opts.bPercID) ){
        Log(&rLog, LOG_FATAL, "Percentage Identity cannot be calculated if Kimura Distances are used.");
    }

    /* iteration destroys effect of pile-up */
    if ( (TRUE == opts->aln_opts.bPileup) && (opts->aln_opts.iNumIterations > 0) ){
        Log(&rLog, LOG_WARN, "Iteration destroys effect of pile-up.");
    }

    AlnOptsLogicCheck(& opts->aln_opts);
}
/* end of UserOptsLogicCheck */



/**
 * @brief Parse command line parameters. Will exit if help/usage etc
 * are called or or call Log(&rLog, LOG_FATAL, ) if an error was detected.
 *
 * @param[out] user_opts
 * User parameter struct, with defaults already set.
 * @param[in] argc
 * mains argc
 * @param[in] argv
 * mains argv
 * 
 */    
void
ParseCommandLine(cmdline_opts_t *user_opts, int argc, char **argv)
{

     /* argtable command line parsing:
     * see
     * http://argtable.sourceforge.net/doc/argtable2-intro.html
     *
     * basic structure is: arg_xxxN:
     * xxx can be int, lit, db1, str, rex, file or date
     * If N = 0, arguments may appear zero-or-once; N = 1 means
     * exactly once, N = n means up to maxcount times
     *
     *
     * @note: changes here, might also affect main.cpp:ConvertOldCmdLine()
     *
     */  
   
    struct arg_rem  *rem_seq_input  = arg_rem(NULL, "\nSequence Input:");
    struct arg_file *opt_seqin = arg_file0("i", "in,infile",
                                            "{<file>,-}",
                                            "Multiple sequence input file (- for stdin)");
    struct arg_file *opt_hmm_in = arg_filen(NULL, "hmm-in", "<file>",
                                            /*min*/ 0, /*max*/ 128,
                                            "HMM input files");
    struct arg_file *opt_hmm_batch = arg_file0(NULL, "hmm-batch", "<file>", /* FS, r291 -> */
                                               "specify HMMs for individual sequences");
    struct arg_lit *opt_dealign = arg_lit0(NULL, "dealign",
                                           "Dealign input sequences");
    struct arg_file *opt_profile1 = arg_file0(NULL, "profile1,p1",
                                              "<file>",
                                              "Pre-aligned multiple sequence file (aligned columns will be kept fix)");
    struct arg_file *opt_profile2 = arg_file0(NULL, "profile2,p2",
                                              "<file>",
                                              "Pre-aligned multiple sequence file (aligned columns will be kept fix)");
    struct arg_lit *opt_isprofile = arg_lit0(NULL, "is-profile",
                                             "disable check if profile, force profile (default no)");

    struct arg_str *opt_seqtype = arg_str0("t", "seqtype",
                                           "{Protein, RNA, DNA}",
                                           "Force a sequence type (default: auto)");
/*    struct arg_lit *opt_force_protein = arg_lit0(NULL, "protein",
                                         "Set sequence type to protein even if Clustal guessed nucleic acid"); */
    struct arg_str *opt_infmt = arg_str0(NULL, "infmt",
                                            "{a2m=fa[sta],clu[stal],msf,phy[lip],selex,st[ockholm],vie[nna]}",
                                            "Forced sequence input file format (default: auto)");
    struct arg_lit *opt_resno = arg_lit0(NULL, "residuenumber,resno",
                                         "in Clustal format print residue numbers (default no)");
    
    struct arg_rem  *rem_guidetree  = arg_rem(NULL, "\nClustering:");
    struct arg_str *opt_pairdist = arg_str0("p", "pairdist",
                                             "{ktuple}",
                                             "Pairwise alignment distance measure");
    struct arg_file *opt_distmat_in = arg_file0(NULL, "distmat-in",
                                                "<file>",
                                                "Pairwise distance matrix input file (skips distance computation)");
    struct arg_file *opt_distmat_out = arg_file0(NULL, "distmat-out",
                                                 "<file>",
                                                 "Pairwise distance matrix output file");
    struct arg_file *opt_guidetree_in = arg_file0(NULL, "guidetree-in",
                                                  "<file>",
                                                  "Guide tree input file (skips distance computation and guide-tree clustering step)");
    struct arg_file *opt_guidetree_out = arg_file0(NULL, "guidetree-out",
                                                   "<file>",
                                                   "Guide tree output file");
    /* AW: mbed is default since at least R253
       struct arg_lit *opt_mbed = arg_lit0(NULL, "mbed",
       "Fast, Mbed-like clustering for guide-tree calculation");
       struct arg_lit *opt_mbed_iter = arg_lit0(NULL, "mbed-iter",
       "Use Mbed-like clustering also during iteration");
    */
    struct arg_lit *opt_pileup = arg_lit0(NULL, "pileup", 
                                          "Sequentially align sequences");
    struct arg_lit *opt_full = arg_lit0(NULL, "full",
                                        "Use full distance matrix for guide-tree calculation (might be slow; mBed is default)");
    struct arg_lit *opt_full_iter = arg_lit0(NULL, "full-iter",
                                        "Use full distance matrix for guide-tree calculation during iteration (might be slowish; mBed is default)");

    struct arg_int *opt_clustersize = arg_int0(NULL, "cluster-size", "<n>", 
                                               "soft maximum of sequences in sub-clusters"); /* FS, r274 -> */

    struct arg_file *opt_clustfile = arg_file0(NULL, "clustering-out",
                                               "<file>",
                                               "Clustering output file"); /* FS, r274 -> */
    struct arg_int *opt_trans = arg_int0(NULL, "trans", "<n>", "use transitivity (default: 0)");
    struct arg_file *opt_posteriorfile = arg_file0(NULL, "posterior-out",
                                               "<file>",
                                               "Posterior probability output file"); /* FS, r288 -> */
    struct arg_lit *opt_usekimura = arg_lit0(NULL, "use-kimura",
                                             "use Kimura distance correction for aligned sequences (default no)");
    struct arg_lit *opt_percentid = arg_lit0(NULL, "percent-id",
                                             "convert distances into percent identities (default no)");

    struct arg_str *opt_clustering = arg_str0("c", "clustering",
                                              "{UPGMA}",
                                              "Clustering method for guide tree");

    
    struct arg_rem *rem_aln_output  = arg_rem(NULL, "\nAlignment Output:");
    struct arg_file *opt_outfile = arg_file0("o", "out,outfile",
                                             "{file,-}",
                                             "Multiple sequence alignment output file (default: stdout)");
    struct arg_str *opt_outfmt = arg_str0(NULL, "outfmt",
                                            "{a2m=fa[sta],clu[stal],msf,phy[lip],selex,st[ockholm],vie[nna]}",
                                            "MSA output file format (default: fasta)");
    struct arg_int *opt_wrap = arg_int0(NULL, "wrap", "<n>", 
                                        "number of residues before line-wrap in output");
    struct arg_str *opt_output_order = arg_str0(NULL, "output-order",
                                                "{input-order,tree-order}",
                                                "MSA output order like in input/guide-tree");

    
    struct arg_rem *rem_iteration  = arg_rem(NULL, "\nIteration:");
    struct arg_str *opt_num_iterations = arg_str0(NULL, "iterations,iter",
                                                  /* FIXME "{<n>,auto}", "Number of combined guide-tree/HMM iterations"); */
                                                  "<n>", "Number of (combined guide-tree/HMM) iterations");
    struct arg_int *opt_max_guidetree_iterations = arg_int0(NULL, "max-guidetree-iterations",
                                                            "<n>", "Maximum number of guidetree iterations");
    struct arg_int *opt_max_hmm_iterations = arg_int0(NULL, "max-hmm-iterations",
                                                      "<n>", "Maximum number of HMM iterations");

   
    struct arg_rem *rem_limits  = arg_rem(NULL, "\nLimits (will exit early, if exceeded):");
    struct arg_int *opt_max_seq = arg_int0(NULL, "maxnumseq", "<n>",
                                           "Maximum allowed number of sequences");
    struct arg_int *opt_max_seqlen = arg_int0(NULL, "maxseqlen", "<l>", 
                                              "Maximum allowed sequence length");


    struct arg_rem *rem_misc  = arg_rem(NULL, "\nMiscellaneous:");

    struct arg_lit *opt_autooptions = arg_lit0(NULL, "auto",
                                         "Set options automatically (might overwrite some of your options)");
    struct arg_int *opt_threads = arg_int0(NULL, "threads", "<n>", 
                                              "Number of processors to use");
    struct arg_file *opt_pseudo = arg_file0(NULL, "pseudo", "<file>",
                                            "Input file for pseudo-count parameters");
    struct arg_file *opt_logfile = arg_file0("l", "log",
                                             "<file>",
                                             "Log all non-essential output to this file");
    struct arg_lit *opt_help = arg_lit0("h", "help",
                                         "Print this help and exit");
    struct arg_lit *opt_version = arg_lit0(NULL, "version",
                                           "Print version information and exit");
    struct arg_lit *opt_long_version = arg_lit0(NULL, "long-version",
                                           "Print long version information and exit");
    struct arg_lit *opt_verbose = arg_litn("v", "verbose",
                                           0, 3,
                                           "Verbose output (increases if given multiple times)");
    struct arg_lit *opt_force = arg_lit0(NULL, "force",
                                         "Force file overwriting");
    struct arg_int *opt_macram = arg_int0(NULL, "MAC-RAM", "<n>", /* keep this quiet for the moment, FS r240 -> */
                                          NULL/*"maximum amount of RAM to use for MAC algorithm (in MB)"*/);

    struct arg_end *opt_end = arg_end(10); /* maximum number of errors
                                            * to store */

    void *argtable[] = {rem_seq_input,
                        opt_seqin,
                        opt_hmm_in,
                        opt_hmm_batch, /* FS, r291 -> */
                        opt_dealign,
                        opt_profile1,
                        opt_profile2,
                        opt_isprofile, /* FS, r282 ->*/
                        opt_seqtype,
                        /* opt_force_protein, */
                        opt_infmt,
                        rem_guidetree,
#if 0
                        /* no other options then default available or not implemented */
                        opt_pairdist,
#endif
                        opt_distmat_in,
                        opt_distmat_out,
                        opt_guidetree_in,
                        opt_guidetree_out,
                        opt_pileup, /* FS, r288 -> */
                        opt_full, /* FS, r250 -> */
                        opt_full_iter, /* FS, r250 -> */
                        opt_clustersize, /* FS, r274 -> */
                        opt_clustfile, /* FS, r274 -> */
                        opt_trans, /* FS, r290 -> */
                        opt_posteriorfile, /* FS, r288 -> */
                        opt_usekimura, /* FS, r282 ->*/
                        opt_percentid, /* FS, r282 ->*/
#if 0
                        /* no other options then default available */
                        opt_clustering,
#endif
                        rem_aln_output,
                        opt_outfile,
                        opt_outfmt,
                        opt_resno,  /* FS, 274 -> */
                        opt_wrap, /* FS, 274 -> */
                        opt_output_order, /* FS, 274 -> */

                        rem_iteration,
                        opt_num_iterations,
                        opt_max_guidetree_iterations,
                        opt_max_hmm_iterations,

                        rem_limits,
                        opt_max_seq,
                        opt_max_seqlen,

                        rem_misc,
                        opt_autooptions,
                        opt_threads,
                        opt_pseudo,
                        opt_logfile,
                        opt_help,
                        opt_verbose,
                        opt_version,
                        opt_long_version,
                        opt_force,
                        opt_macram, /* FS, r240 -> r241 */

                        opt_end};
    int nerrors;


    /* Verify the argtable[] entries were allocated sucessfully
     */
    if (arg_nullcheck(argtable)) {
        Log(&rLog, LOG_FATAL, "insufficient memory (for argtable allocation)");
    }

    /* Parse the command line as defined by argtable[]
     */
    nerrors = arg_parse(argc, argv, argtable);

    /* Special case: '--help' takes precedence over error reporting
     */
    if (opt_help->count > 0) {
        printf("%s - %s (%s)\n", PACKAGE_NAME, PACKAGE_VERSION, PACKAGE_CODENAME);

        printf("\n");
        printf("If you like Clustal-Omega please cite:\n%s\n", CITATION);
        printf("If you don't like Clustal-Omega, please let us know why (and cite us anyway).\n");
        /* printf("You can contact reach us under %s\n", PACKAGE_BUGREPORT); */
        printf("\n");
        printf("Check http://www.clustal.org for more information and updates.\n");
            
        printf("\n");
        printf("Usage: %s", basename(argv[0]));
        arg_print_syntax(stdout,argtable, "\n");

        printf("\n");
        printf("A typical invocation would be: %s -i my-in-seqs.fa -o my-out-seqs.fa -v\n",
               basename(argv[0]));
        printf("See below for a list of all options.\n");

        arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));
        exit(EXIT_SUCCESS);
    }

    /* Special case: '--version' takes precedence over error reporting
     */
    if (opt_version->count > 0) {
        printf("%s\n", PACKAGE_VERSION);
        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
        exit(EXIT_SUCCESS);
    }

    /* Special case: '--long-version' takes precedence over error reporting
     */
    if (opt_long_version->count > 0) {
        char zcLongVersion[1024];
        PrintLongVersion(zcLongVersion, sizeof(zcLongVersion));
        printf("%s\n", zcLongVersion);
        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
        exit(EXIT_SUCCESS);
    }

    /* If the parser returned any errors then display them and exit
     */
    if (nerrors > 0) {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout, opt_end, PACKAGE);
        fprintf(stderr, "For more information try: %s --help\n", argv[0]);
        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
        exit(EXIT_FAILURE);
    }

    
    /* ------------------------------------------------------------
     *
     * Command line successfully parsed. Now transfer values to
     * user_opts. While doing so, make sure that given input files
     * exist and given output files are writable do not exist, or if
     * they do, should be overwritten.
     *
     * No logic checks here! They are done in a different function
     *
     * ------------------------------------------------------------*/
        
    
    /* not part of user_opts because it declared in src/util.h */
    if (0 == opt_verbose->count) {
        rLog.iLogLevelEnabled = LOG_WARN;
    } else if (1 == opt_verbose->count) {
        rLog.iLogLevelEnabled = LOG_INFO;
    } else if (2 == opt_verbose->count) {
        rLog.iLogLevelEnabled = LOG_VERBOSE;
    } else if (3 == opt_verbose->count) {
        rLog.iLogLevelEnabled = LOG_DEBUG;
    }

    user_opts->aln_opts.bAutoOptions = opt_autooptions->count;

    user_opts->bDealignInputSeqs = opt_dealign->count;

    /* NOTE: full distance matrix used to be default - there was
       --mbed flag but no --full flag after r250 decided that mBed
       should be default - now need --full flag to turn off mBed.
       wanted to retain mBed Boolean, so simply added --full flag. if
       both flags set (erroneously) want --mbed to overwrite --full,
       therefore do --full 1st, the --mbed, FS, r250 */
    if (opt_full->count){
        user_opts->aln_opts.bUseMbed = FALSE;
    }

    if (opt_full_iter->count){
        user_opts->aln_opts.bUseMbedForIteration = FALSE;
    }

    if (opt_pileup->count){
        user_opts->aln_opts.bPileup = TRUE;
    }


    user_opts->bForceFileOverwrite = opt_force->count;

    /* log-file
     */
    if (opt_logfile->count > 0) {
        user_opts->pcLogFile = CkStrdup(opt_logfile->filename[0]);
        
        /* warn if already exists or not writable */
        if (FileExists(user_opts->pcLogFile) && ! user_opts->bForceFileOverwrite) {
            Log(&rLog, LOG_FATAL, "%s '%s'. %s",
                  "Cowardly refusing to overwrite already existing file",
                  user_opts->pcLogFile,
                  "Use --force to force overwriting.");
        }
        if (! FileIsWritable(user_opts->pcLogFile)) {
            Log(&rLog, LOG_FATAL, "Sorry, I do not have permission to write to file '%s'.",
                  user_opts->pcLogFile);
        }
    }


    /* pseudo-count-file
     */
    if (opt_pseudo->count > 0) {
        user_opts->pcPseudoFile = CkStrdup(opt_pseudo->filename[0]);
        if (! FileExists(user_opts->pcPseudoFile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->pcPseudoFile);
        }
    }


    /* normal sequence input (no profile)
     */
    if (opt_seqin->count > 0) {
        user_opts->pcSeqInfile = CkStrdup(opt_seqin->filename[0]);
    }

    /* Input limitations
     */
    /* maximum number of sequences */
    if (opt_max_seq->count > 0) {
        user_opts->iMaxNumSeq = opt_max_seq->ival[0];
    }
    
    /* maximum sequence length */
    if (opt_max_seqlen->count > 0) {
        user_opts->iMaxSeqLen = opt_max_seqlen->ival[0];
    }
    
    /* Output format
     */
    if (opt_infmt->count > 0) {
        /* avoid gcc warning about discarded qualifier */
        char *tmp = (char *)opt_infmt->sval[0];
        user_opts->iSeqInFormat = String2SeqfileFormat(tmp);
    } else {
        user_opts->iSeqInFormat = SQFILE_UNKNOWN;
    }


    /* Sequence type
     */
    if (opt_seqtype->count > 0) {
        if (STR_NC_EQ(opt_seqtype->sval[0], "protein")) {
            user_opts->iSeqType = SEQTYPE_PROTEIN;
        } else if (STR_NC_EQ(opt_seqtype->sval[0], "rna")) {
            user_opts->iSeqType = SEQTYPE_RNA;
        } else if  (STR_NC_EQ(opt_seqtype->sval[0], "dna")) {
            user_opts->iSeqType = SEQTYPE_DNA;
        } else {
            Log(&rLog, LOG_FATAL, "Unknown sequence type '%s'", opt_seqtype->sval[0]);
        }
    }
/*    if (opt_force_protein->count > 0) {
        user_opts->iSeqType = SEQTYPE_PROTEIN;
    } */

    /* Profile input
     */
    if (opt_profile1->count > 0) {
        user_opts->pcProfile1Infile = CkStrdup(opt_profile1->filename[0]);
        if (! FileExists(user_opts->pcProfile1Infile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->pcProfile1Infile);
        }
    }
    
    if (opt_profile2->count > 0) {
        user_opts->pcProfile2Infile = CkStrdup(opt_profile2->filename[0]);
        if (! FileExists(user_opts->pcProfile2Infile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->pcProfile2Infile);
        }
    }
    
    if (opt_isprofile->count){
        user_opts->bIsProfile = TRUE; 
    }
    if (opt_usekimura->count){
        user_opts->aln_opts.bUseKimura = TRUE;
    }
    if (opt_percentid->count){
        user_opts->aln_opts.bPercID = TRUE;
    }

    
    /* HMM input
     */
    user_opts->aln_opts.iHMMInputFiles = 0;
    user_opts->aln_opts.ppcHMMInput = NULL;
    if (opt_hmm_in->count>0) {
        int iAux;
        user_opts->aln_opts.iHMMInputFiles = opt_hmm_in->count;
        user_opts->aln_opts.ppcHMMInput = (char **) CKMALLOC(
            user_opts->aln_opts.iHMMInputFiles * sizeof(char*));
        for (iAux=0; iAux<opt_hmm_in->count; iAux++) {
            user_opts->aln_opts.ppcHMMInput[iAux] = CkStrdup(opt_hmm_in->filename[iAux]);
            if (! FileExists(user_opts->aln_opts.ppcHMMInput[iAux])) {
                Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->aln_opts.ppcHMMInput[iAux]);
            }
        }
    }


    /* HMM Batch, FS, r291 ->
     */
    user_opts->aln_opts.pcHMMBatch = NULL;
    if (opt_hmm_batch->count>0){
        user_opts->aln_opts.pcHMMBatch = CkStrdup(opt_hmm_batch->filename[0]);
        if (! FileExists(user_opts->aln_opts.pcHMMBatch)){
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->aln_opts.pcHMMBatch);
        }
    }

    /* Pair distance method
     */
    if (opt_pairdist->count > 0) {
        if (STR_NC_EQ(opt_pairdist->sval[0], "ktuple")) {
            user_opts->aln_opts.iPairDistType = PAIRDIST_KTUPLE;
        } else {
            Log(&rLog, LOG_FATAL, "Unknown pairdist method '%s'", opt_pairdist->sval[0]);
        }
    }


    /* Distance matrix input
     */
    if (opt_distmat_in->count > 0) {
        user_opts->aln_opts.pcDistmatInfile = CkStrdup(opt_distmat_in->filename[0]);
        if (! FileExists(user_opts->aln_opts.pcDistmatInfile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->aln_opts.pcDistmatInfile);
        }
    }


    /* Distance matrix output
     */
    if (opt_distmat_out->count > 0) {
        user_opts->aln_opts.pcDistmatOutfile = CkStrdup(opt_distmat_out->filename[0]);
        
        /* warn if already exists or not writable */
        if (FileExists(user_opts->aln_opts.pcDistmatOutfile) && ! user_opts->bForceFileOverwrite) {
            Log(&rLog, LOG_FATAL, "%s '%s'. %s",
                  "Cowardly refusing to overwrite already existing file",
                  user_opts->aln_opts.pcDistmatOutfile,
                  "Use --force to force overwriting.");
        }
        if (! FileIsWritable(user_opts->aln_opts.pcDistmatOutfile)) {
            Log(&rLog, LOG_FATAL, "Sorry, I do not have permission to write to file '%s'.",
                user_opts->aln_opts.pcDistmatOutfile);
        }
    }

    /* Clustering
     *
     */
    if (opt_clustering->count > 0) {
        if (STR_NC_EQ(opt_clustering->sval[0], "upgma")) {
            user_opts->aln_opts.iClusteringType = CLUSTERING_UPGMA;
        } else {
            Log(&rLog, LOG_FATAL, "Unknown guide-tree clustering method '%s'", opt_clustering->sval[0]);
        }
    }

    if (opt_clustersize->count > 0){ /* FS, r274 -> */
        if (opt_clustersize->ival[0] > 0){
            user_opts->aln_opts.iClustersizes = opt_clustersize->ival[0];
        }
    }

    if (opt_clustfile->count > 0){ /* FS, r274 -> */
        user_opts->aln_opts.pcClustfile = CkStrdup(opt_clustfile->filename[0]);
        /*if (! FileExists(user_opts->aln_opts.pcClustfile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->aln_opts.pcClustfile);
            }*/
    }

    if (opt_trans->count > 0){ /* FS, r290 -> */
        user_opts->aln_opts.iTransitivity = opt_trans->ival[0];
    }
    if (opt_posteriorfile->count > 0){ /* FS, r288 -> */
        user_opts->aln_opts.pcPosteriorfile = CkStrdup(opt_posteriorfile->filename[0]);
        /*if (! FileExists(user_opts->aln_opts.pcPosteriorfile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->aln_opts.pcPosteriorfile);
            }*/
    }


    /* Guidetree input
     */
    if (opt_guidetree_in->count > 0) {
        user_opts->aln_opts.pcGuidetreeInfile = CkStrdup(opt_guidetree_in->filename[0]);
        if (! FileExists(user_opts->aln_opts.pcGuidetreeInfile)) {
            Log(&rLog, LOG_FATAL, "File '%s' does not exist.", user_opts->aln_opts.pcGuidetreeInfile);
        }
    }
    
    
    /* Guidetree output
     */
    if (opt_guidetree_out->count > 0) {
        user_opts->aln_opts.pcGuidetreeOutfile = CkStrdup(opt_guidetree_out->filename[0]);
        
        /* warn if already exists or not writable */
        if (FileExists(user_opts->aln_opts.pcGuidetreeOutfile) && ! user_opts->bForceFileOverwrite) {
            Log(&rLog, LOG_FATAL, "%s '%s'. %s",
                  "Cowardly refusing to overwrite already existing file",
                  user_opts->aln_opts.pcGuidetreeOutfile,
                  "Use --force to force overwriting.");
        }
        if (! FileIsWritable(user_opts->aln_opts.pcGuidetreeOutfile)) {
            Log(&rLog, LOG_FATAL, "Sorry, I do not have permission to write to file '%s'.",
                  user_opts->aln_opts.pcGuidetreeOutfile);
        }
    }
    

    /* max guidetree iterations
     */
    if (opt_max_guidetree_iterations->count > 0) {
        user_opts->aln_opts.iMaxGuidetreeIterations = opt_max_guidetree_iterations->ival[0];
    }


    /* max guidetree iterations
     */
    if (opt_max_hmm_iterations->count > 0) {
        user_opts->aln_opts.iMaxHMMIterations = opt_max_hmm_iterations->ival[0];
    }

     /* number of iterations
      */
     if (opt_num_iterations->count > 0) {
        if (STR_NC_EQ(opt_num_iterations->sval[0], "auto")) {
            Log(&rLog, LOG_FATAL, "Automatic iteration not supported at the moment.");
            user_opts->aln_opts.bIterationsAuto = TRUE;

        } else {
            int iAux;
            user_opts->aln_opts.bIterationsAuto = FALSE;
            for (iAux=0; iAux<(int)strlen(opt_num_iterations->sval[0]); iAux++) {
                if (! isdigit(opt_num_iterations->sval[0][iAux])) {
                    Log(&rLog, LOG_FATAL, "Couldn't iteration parameter: %s",
                          opt_num_iterations->sval[0]);
                }
            }
            user_opts->aln_opts.iNumIterations = atoi(opt_num_iterations->sval[0]);
        }
    }

    
    /* Alignment output
     */
    if (opt_outfile->count > 0) {
        user_opts->pcAlnOutfile = CkStrdup(opt_outfile->filename[0]);

        /* warn if already exists or not writable */
        if (FileExists(user_opts->pcAlnOutfile) && ! user_opts->bForceFileOverwrite) {
            Log(&rLog, LOG_FATAL, "%s '%s'. %s",
                  "Cowardly refusing to overwrite already existing file",
                  user_opts->pcAlnOutfile,
                  "Use --force to force overwriting.");
        }
        if (! FileIsWritable(user_opts->pcAlnOutfile)) {
            Log(&rLog, LOG_FATAL, "Sorry, I do not have permission to write to file '%s'.",
                  user_opts->pcAlnOutfile);
        }
    }
    

    /* Output format
     */
    if (opt_outfmt->count > 0) {
        /* avoid gcc warning about discarded qualifier */
        char *tmp = (char *)opt_outfmt->sval[0];
        user_opts->iAlnOutFormat = String2SeqfileFormat(tmp);
        if (SQFILE_UNKNOWN == user_opts->iAlnOutFormat) {
            Log(&rLog, LOG_FATAL, "Unknown output format '%s'", opt_outfmt->sval[0]);
        }
    }

    /* Number of threads
     */
#ifdef HAVE_OPENMP
    if (opt_threads->count > 0) {
        if (opt_threads->ival[0] <= 0) {
            Log(&rLog, LOG_FATAL, "Changing number of threads to %d doesn't make sense.", 
                  opt_threads->ival[0]);    
        }
        user_opts->iThreads = opt_threads->ival[0];
    }

#else
    if (opt_threads->count > 0) {
        if (opt_threads->ival[0] > 1) {
            Log(&rLog, LOG_FATAL, "Cannot change number of threads to %d. %s was build without OpenMP support.", 
                  opt_threads->ival[0], PACKAGE_NAME);    
        }
    }
#endif


    /* max MAC RAM (maximum amount of RAM set aside for MAC algorithm)
     */
    if (opt_macram->count > 0) { /* FS, r240 -> r241 */
        user_opts->aln_opts.rHhalignPara.iMacRamMB = opt_macram->ival[0];
    }

    /* Number of residues in output before line-wrap
     */
    if (opt_wrap->count > 0) { /* FS, r274 -> */
        user_opts->iWrap = opt_wrap->ival[0];
    }

    user_opts->bResno = opt_resno->count;

    /* output-order
     * like input (INPUT_ORDER = 0) or tree (TREE_ORDER = 1)
     * if output-order format not valid use INPUT_ORDER
     */
    if (opt_output_order->count > 0){
        user_opts->iOutputOrder = (0 == strcmp(opt_output_order->sval[0], "input-order")) ? INPUT_ORDER : 
            (0 == strcmp(opt_output_order->sval[0], "tree-order")) ? TREE_ORDER : INPUT_ORDER;
    }


    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    UserOptsLogicCheck(user_opts);

    return; 
}
/* end of ParseCommandLine() */




/**
 *
 * @brief the 'real' main function
 *
 */
int
MyMain(int argc, char **argv)
{
    mseq_t *prMSeq = NULL;
    mseq_t *prMSeqProfile1 = NULL;
    mseq_t *prMSeqProfile2 = NULL;
    cmdline_opts_t cmdline_opts = {0};

    /* Must happen first: setup logger */
    LogDefaultSetup(&rLog);

    /*Log(&rLog, LOG_WARN, "This is a non-public realase of %s. Please do not distribute.", PACKAGE_NAME);*/
    /*Log(&rLog, LOG_WARN, "This is a beta version of %s, for protein only.", PACKAGE_NAME);*/ /* FS, r237 -> 238 */

    SetDefaultUserOpts(&(cmdline_opts));

    ParseCommandLine(&cmdline_opts, argc, argv);
    
    if (NULL != cmdline_opts.pcLogFile) {
        prLogFile = fopen(cmdline_opts.pcLogFile, "w");
        LogSetFP(&rLog, LOG_INFO, prLogFile);
        LogSetFP(&rLog, LOG_VERBOSE, prLogFile);
        LogSetFP(&rLog, LOG_DEBUG, prLogFile);
    }

    InitClustalOmega(cmdline_opts.iThreads);

    if (rLog.iLogLevelEnabled < LOG_INFO) {
        PrintUserOpts(LogGetFP(&rLog, LOG_INFO), & cmdline_opts);
        PrintAlnOpts(LogGetFP(&rLog, LOG_INFO), & (cmdline_opts.aln_opts));
    }
#if 1
    if (NULL != cmdline_opts.pcPseudoFile){
        ReadPseudoCountParams(&cmdline_opts.aln_opts.rHhalignPara, cmdline_opts.pcPseudoFile);
    }
#endif

    /* Read sequence file
     *
     */
    if (NULL != cmdline_opts.pcSeqInfile) {
        NewMSeq(&prMSeq);
        if (ReadSequences(prMSeq, cmdline_opts.pcSeqInfile,
                          cmdline_opts.iSeqType, cmdline_opts.iSeqInFormat,
                          cmdline_opts.bIsProfile, cmdline_opts.bDealignInputSeqs, 
                          cmdline_opts.iMaxNumSeq, cmdline_opts.iMaxSeqLen, cmdline_opts.aln_opts.pcHMMBatch)) {
            Log(&rLog, LOG_FATAL, "Reading sequence file '%s' failed", cmdline_opts.pcSeqInfile);
        }
#if TRACE
        {
            int iAux;
            for (iAux=0; iAux<prMSeq->nseqs; iAux++) {
                Log(&rLog, LOG_FORCED_DEBUG, "seq no %d: seq = %s", iAux, prMSeq->seq[iAux]);
                LogSqInfo(&prMSeq->sqinfo[iAux]);
            }
        }
#endif
    }
    /* k-tuple pairwise distance calculation seg-faults if 
     * only one sequence, simply exit early.
     * note that for profile/profile alignment prMSeq is NULL 
     * FS, r222->r223 */
    if (prMSeq && (prMSeq->nseqs <= 1)){
        Log(&rLog, LOG_FATAL, "File '%s' contains %d sequence%s, nothing to align",
              cmdline_opts.pcSeqInfile, prMSeq->nseqs, 1==prMSeq->nseqs?"":"s");
    }
    /* if there are fewer sequences than target size of clusters,
     * then mBed is unnecessary, FS, r283-> */
    if ( (prMSeq) && (prMSeq->nseqs <= cmdline_opts.aln_opts.iClustersizes) ){
        cmdline_opts.aln_opts.bUseMbed = FALSE;
        cmdline_opts.aln_opts.bUseMbedForIteration = FALSE;
        Log(&rLog, LOG_INFO, "not more sequences (%d) than cluster-size (%d), turn off mBed", 
            prMSeq->nseqs, cmdline_opts.aln_opts.iClustersizes);
    }

    /* Dealign if requested and neccessary
     */
    if (NULL != prMSeq) {
        if (/*TRUE == prMSeq->aligned &&*/ cmdline_opts.bDealignInputSeqs) {
            Log(&rLog, LOG_INFO, "Dealigning already aligned input sequences as requested.");
            DealignMSeq(prMSeq);
        }
    }


    /* Read profile1
     *
     */
    if (NULL != cmdline_opts.pcProfile1Infile) {
        NewMSeq(&prMSeqProfile1);
        if (ReadSequences(prMSeqProfile1, cmdline_opts.pcProfile1Infile,
                          cmdline_opts.iSeqType,  cmdline_opts.iSeqInFormat,
                          cmdline_opts.bIsProfile, cmdline_opts.bDealignInputSeqs, 
                          cmdline_opts.iMaxNumSeq, cmdline_opts.iMaxSeqLen, cmdline_opts.aln_opts.pcHMMBatch)) {
            Log(&rLog, LOG_FATAL, "Reading sequences from profile file '%s' failed",
                  cmdline_opts.pcProfile1Infile);
        }
        /* FIXME: commented out. FS, r240 -> r241  
         * for explanation see below */
        /*if (1==prMSeqProfile1->nseqs) {
          Log(&rLog, LOG_FATAL, "'%s' contains only one sequence and can therefore not be used as a profile",
          cmdline_opts.pcProfile1Infile);
          }*/
        if (FALSE == prMSeqProfile1->aligned) {
            Log(&rLog, LOG_FATAL, "Sequences in '%s' are not aligned, i.e. this is not a profile",
                  cmdline_opts.pcProfile1Infile);
        }
    }

    

    /* Read profile2
     *
     */
    if (NULL != cmdline_opts.pcProfile2Infile) {
        NewMSeq(&prMSeqProfile2);
        if (ReadSequences(prMSeqProfile2, cmdline_opts.pcProfile2Infile,
                          cmdline_opts.iSeqType,  cmdline_opts.iSeqInFormat,
                          cmdline_opts.bIsProfile, cmdline_opts.bDealignInputSeqs, 
                          cmdline_opts.iMaxNumSeq, cmdline_opts.iMaxSeqLen, cmdline_opts.aln_opts.pcHMMBatch)) {
            Log(&rLog, LOG_FATAL, "Reading sequences from profile file '%s' failed",
                  cmdline_opts.pcProfile2Infile);
        }
        /* FIXME: there is no (clean) way to align a single sequence to a profile. 
         * if we go down the -i route, it causes a seg-fault in the pair-wise 
         * k-tuple distance calculation. However, single sequences can be 
         * understood as 1-profiles. Therefore we have to allow for 1-profiles.
         * FS, r240 -> r241 
         */
        /*if (1==prMSeqProfile2->nseqs) {
          Log(&rLog, LOG_FATAL, "'%s' contains only one sequence and can therefore not be used as a profile",
          cmdline_opts.pcProfile2Infile);
          }*/
        if (FALSE == prMSeqProfile2->aligned) {
            Log(&rLog, LOG_FATAL, "Sequences in '%s' are not aligned, i.e. this is not a profile",
                  cmdline_opts.pcProfile2Infile);
        }
    }


    /* Depending on the input we got perform
     *
     * (i) normal alignment: seq + optional profile
     * or
     * (ii) profile profile alignment
     *
     */

    if (NULL != prMSeq) {

        if (2 == prMSeq->nseqs){
            /* if there are only 2 sequences then the order does not matter */
            /* this is important, because for pair-wise alignment we don't do 
             * tree indexing*, and trying to use tree-indexing during output 
             * will cause a segmentation fault */
            cmdline_opts.iOutputOrder = INPUT_ORDER;
        }
        if (TREE_ORDER == cmdline_opts.iOutputOrder){
            /* this is a crutch, if tree_order==NULL use input order, 
             * otherwise use guide-tree-order */
            prMSeq->tree_order = (int *)CKMALLOC(prMSeq->nseqs * sizeof(int));
        }
        if (Align(prMSeq, prMSeqProfile1, & cmdline_opts.aln_opts)) {
            Log(&rLog, LOG_FATAL, "An error occured during the alignment");
        }

        if (cmdline_opts.aln_opts.iMaxHMMIterations >= 0){
            if (WriteAlignment(prMSeq, cmdline_opts.pcAlnOutfile, 
                               cmdline_opts.iAlnOutFormat, cmdline_opts.iWrap, cmdline_opts.bResno)) {
                Log(&rLog, LOG_FATAL, "Could not save alignment to %s", cmdline_opts.pcAlnOutfile);
            }
        }
#if 0
        {
            bool bSampling = FALSE; /* better set to TRUE for many sequences */
            bool bReportAll = TRUE;
            AliStat(prMSeq, bSampling, bReportAll);
        }
#endif
        

    } else if (NULL != prMSeqProfile1 && NULL != prMSeqProfile2) {
        if (AlignProfiles(prMSeqProfile1, prMSeqProfile2, 
                          cmdline_opts.aln_opts.rHhalignPara)) {
            Log(&rLog, LOG_FATAL, "An error occured during the alignment");
        }
        if (WriteAlignment(prMSeqProfile1, cmdline_opts.pcAlnOutfile, 
                           cmdline_opts.iAlnOutFormat, cmdline_opts.iWrap, cmdline_opts.bResno)) {
            Log(&rLog, LOG_FATAL, "Could not save alignment to %s", cmdline_opts.pcAlnOutfile);
        }
    }


    /* cleanup
     */
    if (NULL != prMSeq) {
        FreeMSeq(&prMSeq);
    }
    if (NULL != prMSeqProfile1) {
        FreeMSeq(&prMSeqProfile1);
    }
    if (NULL != prMSeqProfile2) {
        FreeMSeq(&prMSeqProfile2);
    }

    FreeUserOpts(&cmdline_opts);

    Log(&rLog, LOG_DEBUG, "Successful program exit");

    if (NULL != cmdline_opts.pcLogFile) {
        fclose(prLogFile);
    }
    return EXIT_SUCCESS;
}
/*  end of MyMain() */


