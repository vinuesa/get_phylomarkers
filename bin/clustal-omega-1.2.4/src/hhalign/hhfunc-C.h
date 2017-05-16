/* -*- mode: c; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

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
 *  RCS $Id: hhfunc-C.h 309 2016-06-13 14:10:11Z fabian $
 */

/*
 * Changelog: Michael Remmert made changes to hhalign stand-alone code 
 * FS implemented some of the changes on 2010-11-11 -> MR1
 *
 * did not incorporate SS prediction PSIpred (yet); functions are: 
 * CalculateSS(3), 
 */


// hhfunc.C

//#include "new_new.h" /* memory tracking */

/**
 * AddBackgroundFrequencies()
 *
 * @brief add background frequencies (derived from HMM) to 
 * sequence/profile
 *
 * @param[in,out] ppfFreq,
 * [in] residue frequencies of sequence/profile, 
 * [out] overlayed with HMM background frequencies
 * @param[in,out] ppfPseudoF,
 * [in] residue frequencies+pseudocounts of sequence/profile, 
 * [out] overlayed with HMM background frequencies+pseudocounts
 * @param[in] iSeqLen,
 * length of sequence/profile (not aligned to HMM)
 * @pram[in] prHMM, 
 * background HMM
 * @param[in] ppcSeq,
 * sequences/profile to be 'softened'
 * @param[in] pcPrealigned,
 * sequence aligned to HMM, this is not quite a consensus,
 * it is identical to 1st sequence but over-writes gap information, 
 * if other sequences in profile have (non-gap) residues
 * @param[in] iPreCnt
 * number of sequences pre-aligned (pcPrealigned is 'consensus' of these sequences) 
 * @param[in] pcRepresent
 * sequence representative of HMM, aligned to pcSeq0
 */
void 
AddBackgroundFrequencies(float **ppfFreq, float **ppfPseudoF, float **ppfPseudoTr, 
                         int iSeqLen, hmm_light *prHMM,
                         char **ppcSeq, char *pcPrealigned, int iPreCnt, char *pcRepresent)
{

    char *pcS = pcPrealigned; /* current residue in pre-aligned sequence */
    char *pcH = pcRepresent;  /* current residue in pre-aligned HMM */
    int iS = 0; /* position in un-aligned sequence (corresponds to pcS) */
    int iH = 0; /* position in un-aligned HMM (corresponds to pcH) */
    int iA; /* residue iterator */
    //int iT; /* transition state iterator */
    float fFWeight = 0.50 / sqrt((float)(iPreCnt)); /* weight of backgroud frequencies */ /* FIXME: tune value, 0.50 default */
    //float fFWeight = 0.75;
    float fFThgiew = UNITY - fFWeight; /* weight of 'true' frequencies */
    //float fGWeight = 0.50 / sqrt((float)(iPreCnt)); /* weight of backgroud transitions */ /* FIXME: tune value, 0.50 default */
    //float fGWeight = 0.50 /*/ (float)(iPreCnt)*/; /* weight of backgroud transitions */ /* FIXME: tune value, 0.50 default */
    //float fGThgiew = UNITY - fGWeight; /* weight of 'true' transitions */
    float fAux; 

    if ( (NULL == pcPrealigned) || (NULL == pcRepresent) ){
        /*printf("%s/%s:%d: WARNING HMM=NULL -- didn't think I would get here (carry on, no danger)\n",
          __FUNCTION__, __FILE__, __LINE__);*/
        return;
    }

    if (NULL == prHMM->p){
        printf("%s:%s:%d: WARNING -- Background Pseudocounts point to NULL\n"
               "\tthis is not intended - don't add background but CONTINUE\n", 
               __FUNCTION__, __FILE__, __LINE__);
        return;

    }
    
    /* FIXME: should be 0 (FS thinks) but -1 gives better results */
    iH = iS = 0/*-1*//*+1*/; 
    while ( ('\0' != *pcS) && ('\0' != *pcH) ){
        
        if ( ('-' != *pcH) && ('-' != *pcS) ){
            iH++;
            iS++;
     
#if 0       
            /* match state 
             * - HMM had a gap in previous position (now closed)
             * FIXME: this does not really work */
            if ((iH > 0) && ('-' == *(pcH-1))){

                for (iT = 0; iT < STATE_TRANSITIONS; iT++){
                    ppfPseudoTr[iS-1][iT] = log2f(zf1SeqRevrt[iT]);
                }
            }
#endif

#if 1
            /* do frequencies -- this is not really useful; 
               frequencies are derived from HMM 
               and contain already pseudocounts (PCs), 
               adding frequencies and then PCs will add PCs _twice_
               results are better than not to add them, 
               but not as good as PCs */
            for (iA = 0; iA < AMINOACIDS; iA++){
                fAux = fFThgiew * ppfFreq[iS][iA] + fFWeight * prHMM->f[iH][iA];
                ppfFreq[iS][iA] = fAux;
            }
#endif
            /* do pseudo-counts */
            for (iA = 0; iA < AMINOACIDS; iA++){
                fAux = 
                    fFThgiew * ppfPseudoF[iS][iA] + fFWeight * prHMM->p[iH][iA];
                ppfPseudoF[iS][iA] = fAux;
            }
#if 0 /* do state transitions */
            for (iT = 0; iT < STATE_TRANSITIONS; iT++){
#if 1
                /* this is a very crude method, 
                   which averages the logarithms of the transitions, 
                   effectively performing a geometric mean - 
                   this presumably violates normalisation.
                   however, results are surprisingly good */
                fAux = 
                    fGThgiew * ppfPseudoTr[iS][iT] + 
                    fGWeight * prHMM->tr[iH][iT];
                ppfPseudoTr[iS][iT] = fAux;
#else /* crude averaging */
                /* There are 2 things to consider:
                   (1) one should really blend the probabilities of 
                   the transitions, however, by default we have 
                   the logarithms thereof. 
                   So must exponentiate, blend, and then take log again.
                   This is expensive, and does not seem to lead to better 
                   results than blending the logarithms 
                   (and violating normalisation)
                   (2) transition probabilities for a single sequence 
                   are really easy, there are no insert/delete transitions. 
                   However, there is a begin state that is different 
                   from the main body. 
                   But again, this does not seem to make a blind bit 
                   of difference
                 */
                if (iS > 0){
                    fAux = 
                        fGThgiew * zf1SeqTrans[iT] + 
                        fGWeight * prHMM->linTr[iH][iT];
                    ppfPseudoTr[iS][iT] = log2f(fAux);
                }
                else {
                    fAux = 
                        fGThgiew * zf1SeqInit[iT]  + 
                        fGWeight * prHMM->linTr[iH][iT];
                    ppfPseudoTr[iS][iT] = log2f(fAux);
                }
#endif /* mixing of linear/log */
            } /* 0 <= iT < STATE_TRANSITIONS */
#endif /* did state transitions */
        } /* was Match -- neither HMM nor sequence have gap */
        
        else if ('-' == *pcH){
            /* sequence opened up gap in HMM */
#if 0
            if ((iH > 0) && ('-' != *(pcH-1)) && (iS > 0)){
                /* this is the first gap in HMM
                 * FIXME: this does not really work */
                for (iT = 0; iT < STATE_TRANSITIONS; iT++){
                    ppfPseudoTr[iS-1][iT] = log2f(zf1SeqDel[iT]);
                }
            }
            else {
                /* do nothing, keep single sequence values exactly as they are*/
            }
#endif
            iS++;
        }
        else if ('-' == *pcS){
            /* here the single sequence has a gap, 
               and the HMM (as a whole) does not. There may be individual gaps 
               in the HMM at this stage. By ignoring this we say that the HMM
               dominates the overall behaviour - as in the M2M state as well
             */
            iH++;
        }

        pcH++; 
        pcS++; 
        
    } /* !EO[seq/hmm] */
    
    return;
    
} /* this is the end of AddBackgroundFrequencies() */



/**
 * ReadAndPrepare()
 *
 * @brief Read HMM input file or transfer alignment 
 * and add pseudocounts etc.
 *
 * @param[in] iRnPtype
 * type of read/prepare 
 * enum {INTERN_ALN_2_HMM = 0, READ_ALN_2_HMM, READ_HMM_2_HMM, INTERN_HMM_2_HMM};
 * @param[in] ppcProf 
 * alignment
 * @param[in] iCnt
 * number of seqs in alignment
 * @param[in,out] prHMM,
 * [in] if sequences read/prepared, [out] if HMM from file
 * @param[in] pcPrealigned,
 * (single) sequence aligned to background HMM
 * @param[in] pcRepresent,
 * sequence representing HMM aligned to individual sequence
 * param[in] pdExWeights
 * (external) sequence weights, derived from tree
 * @param[in] infile
 * name of file with HMM info (formerly also alignment)
 * @param[out] q
 * HMM structure with transition probabilities, residue frequencies etc
 * @param[???] qali
 * FIXME: what is qali?
 *
 * @return FAILURE on error, OK otherwise
 */
int
ReadAndPrepare(int iRnPtype, 
               char **ppcProf, int iCnt, hmm_light *prHMM, 
               char *pcPrealigned, char *pcRepresent, double *pdExWeights, 
               char* infile, HMM& q, Alignment* qali=NULL) {
  
  //#ifndef CLUSTALO_NOFILE
  char path[NAMELEN];
  
  /* NOTE: there are different scenarios:

     (i) ("" != infile) - read HMM from file, 
     transfer frequency/transition/pseudo-count (f/t/p) info into prHMM

     (ii) ('\0' != ppcProf[0]) - transfer sequence/alignment into q/qali,
     don't save f/t/p into prHMM, 
     on the contrary, if prior f/t/p available then add it to q/qali, 
     this is only done if (1==iCnt)

     (iii) ('\0' == ppcProf[0]) - re-cycle old HMM information
     recreate a HMM from previous data
  */

  /********************************/
  /*** (o) Recycle internal HMM ***/
  /********************************/
  if ( (INTERN_ALN_2_HMM == iRnPtype) && (iCnt <= 0) ){

      /* NOTE: here we are writing into a HMM structure/class;
         memory has been allocated for this in hhalign.cpp;
         however, as iCnt<=0, there may not be memory for 
         prHMM->n_display sequences/names. 
         But then, there doesn't have to be. 
         At this stage we are just copying one HMM into another HMM, 
         sequences are irrelevant. The only sequence of (marginal) 
         interest is the consensus sequence */
      /* FIXME: check that prHMM is valid -- how? */

      const int ciCons = 0;
      const int ciNoof = ciCons+1;
      const int ciInvalid = -1;
      q.n_display = ciNoof; /* only consensus */
      q.sname = NULL;
      q.ncons      = ciCons;  
      q.nfirst     = ciCons;//prHMM->nfirst;     
      q.nss_dssp   = ciInvalid;//prHMM->nss_dssp;   
      q.nsa_dssp   = ciInvalid;//prHMM->nsa_dssp;   
      q.nss_pred   = ciInvalid;//prHMM->nss_pred;   
      q.nss_conf   = ciInvalid;//prHMM->nss_conf;   
      q.L          = prHMM->L;          
      q.N_in       = prHMM->N_in;       
      q.N_filtered = prHMM->N_filtered; 
      /* NOTE: I (FS) think that only ever 1 sequence will be transferred here,
         that is, the consensus sequence. However, we might want to allow 
         (in the future) to transfer more sequences, 
         hence the awkward for() loop */
#if 0
      for (int i = prHMM->ncons+0; i < prHMM->ncons+q.n_display; i++){
          /* NOTE: In the original hhalign code the first position
             is kept open ('\0'). This makes it difficult to use  
             string functions like strlen/strdup etc.
             Insert a temporary gap (.) to facilitate string operations */
          char cHead  = prHMM->seq[i][0];
          prHMM->seq[i][0] = '.';
          q.seq[i] = strdup(prHMM->seq[i]);
          prHMM->seq[i][0] = q.seq[i][0] = cHead;
      }
#else
      {
          char cHead  = prHMM->seq[prHMM->ncons][0];
          prHMM->seq[prHMM->ncons][0] = '.';
          q.seq[q.ncons] = strdup(prHMM->seq[prHMM->ncons]);
          prHMM->seq[prHMM->ncons][0] = q.seq[q.ncons][0] = cHead;
      }
#endif
      for (int i = 0; i < prHMM->L+1; i++){
          q.Neff_M[i] = prHMM->Neff_M[i];
          q.Neff_I[i] = prHMM->Neff_I[i];
          q.Neff_D[i] = prHMM->Neff_D[i];
      }
      q.Neff_HMM = prHMM->Neff_HMM;
      /* skip longname,name,file,fam,sfam,fold,cl */
      q.lamda = prHMM->lamda;
      q.mu    = prHMM->mu;

      HMMshadow rShadow = {0}; /* friend of HMM to access private members */
      rShadow.copyShadowToHMM(q, *prHMM);

      /* skip trans_lin,ss_dssp,sa_dssp,ss_pred,ss_conf,Xcons */
      /* pav already done in copyShadowToHMM */
      /* skip pnul */

      return OK;


  } /* INTERN_ALN_2_HMM && iCnt<=0 */

  /******************************/
  /*** (i) Read HMM from file ***/
  /******************************/
  char line[LINELEN] = {0}; // input line
  FILE *inf = NULL;  
  //if ( (0 != strcmp(infile,"")) /*&& (iCnt > 0)*/ ) 
  if ( (READ_HMM_2_HMM == iRnPtype) || (READ_ALN_2_HMM == iRnPtype)) {

      if (0 == strcmp(infile,"")){
          printf("%s:%s:%d:\n"
                 "\texpected to re %s from file but no file specified\n"
                 ""
                 , __FUNCTION__, __FILE__, __LINE__
                 , (READ_HMM_2_HMM==iRnPtype?"HMM":"alignment"));
          return FAILURE;
      }
    inf = fopen(infile, "r");
    if (!inf) OpenFileError(infile);
    Pathname(path,infile);

    //}
    //else {
    //inf = stdin;
    //if (v>=2) printf("Reading HMM / multiple alignment from standard input ...\n(To get a help list instead, quit and type %s -h.)\n",program_name);
    //*path='\0';
    //} 

    fgetline(line,LINELEN-1,inf);  
  }

  //if  ( (0 != strcmp(infile,"")) && (iCnt > 0) )
  if ( (READ_HMM_2_HMM == iRnPtype) ){
      
      if (0 == strcmp(infile,"")){
          printf("%s:%s:%d: expected to read HMM from file but no file-name\n", 
                 __FUNCTION__, __FILE__, __LINE__);
      }
      
      // Is infile a HMMER3 file? /* MR1 */
      if (!strncmp(line,"HMMER3",6))
          {
              if (v>=2) {
                  cout<<"Query file is in HMMER3 format\n";
                  cout<<"WARNING! Use of HMMER3 format as input results in dramatically loss of sensitivity!\n";
              }
              
              // Read 'query HMMER file
              rewind(inf);
              q.ReadHMMer3(inf,path);

              // Don't add transition pseudocounts to query!!
              
              // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a] 
              q.PreparePseudocounts();
              
              // DON'T ADD amino acid pseudocounts to query: pcm=0!  q.p[i][a] = f[i][a]
              q.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
              
              /* further down there is a q.Log2LinTransitionProbs 
                 but only if (iCnt>0), however, we still need it it here (i think), 
                 there is no danger of doing this twice, as trans_lin is checked 
                 FIXME (FS, 2011-01-12) */
              /* further down there is a q.Log2LinTransitionProbs 
                 but only if (iCnt>0), however, we still need it it here (i think), 
                 there is no danger of doing this twice, as trans_lin is checked 
                 FIXME (FS, 2011-01-12) */
              //if (par.forward >= 1) 
              {
                  q.Log2LinTransitionProbs(1.0);
              }

          }
      // ... or Is infile an old HMMER file?
      else if (!strncmp(line,"HMMER",5)) {
          if (v>=2) {
              cout<<"Query file is in HMMER format\n";
              cout<<"WARNING! Use of HMMER format as input results in dramatically loss of sensitivity!\n";
          }
          
          // Read 'query HMMER file
          q.ReadHMMer(inf,path);
          if (v>=2 && q.Neff_HMM>11.0) 
              fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);
          
          // Don't add transition pseudocounts to query!!
          
          // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a] 
          q.PreparePseudocounts();
          
          // DON'T ADD amino acid pseudocounts to query: pcm=0!  q.p[i][a] = f[i][a]
          q.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);

          /* further down there is a q.Log2LinTransitionProbs 
             but only if (iCnt>0), however, we still need it it here (i think), 
             there is no danger of doing this twice, as trans_lin is checked 
             FIXME (FS, 2011-01-12) */
          //if (par.forward >= 1) 
              {
                  q.Log2LinTransitionProbs(1.0);
              }

      } /* it was a HMMer file */
      
      // ... or is it an hhm file?
      else if (!strncmp(line,"NAME",4) || !strncmp(line,"HH",2)) {
          
          if (v>=2) cout<<"Query file is in HHM format\n";
          
          // Rewind to beginning of line and read query hhm file
          rewind(inf);
          q.Read(inf,path);
          if (v>=2 && q.Neff_HMM>11.0) 
              fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);
          
          // Add transition pseudocounts to query -> q.p[i][a]
          q.AddTransitionPseudocounts();
          
          // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a] 
          q.PreparePseudocounts();
          
          // Add amino acid pseudocounts to query:  q.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
          q.AddAminoAcidPseudocounts();
          
      } /* it was a HHM file */
      else {
          fprintf(stderr, "%s:%s:%d: Unknown HMM format in infile\n"
                  "infile=%s, #seq=%d\n"
                  , __FUNCTION__, __FILE__, __LINE__
                  , infile, iCnt);
          return FAILURE;
      }
      
      /*fclose(inf);*/
      
      /*** transfer class info to struct */
      prHMM->n_display = q.n_display;
      /* ignore sname*/
      prHMM->seq = (char **)calloc((q.n_display+1), sizeof(char *));
      /* FIXME valgrind says bytes get lost in the above calloc during
       * hmm-iteration
       */
      for (int i = 0; i < q.n_display; i++){
          /* NOTE: In the original hhalign code the first position 
             is kept open ('\0'). This makes it difficult to use 
             string functions like strlen/strdup etc. 
             Insert a temporary gap (.) to facilitate string operations */
          char cHead  = q.seq[i][0];
          q.seq[i][0] = '.';
          prHMM->seq[i] = strdup(q.seq[i]);
          q.seq[i][0] = prHMM->seq[i][0] = cHead; 
      }
      prHMM->ncons      = q.ncons;
      prHMM->nfirst     = q.nfirst;
      prHMM->nss_dssp   = q.nss_dssp;
      prHMM->nsa_dssp   = q.nsa_dssp;
      prHMM->nss_pred   = q.nss_pred;
      prHMM->nss_conf   = q.nss_conf;
      prHMM->L          = q.L;
      prHMM->N_in       = q.N_in;
      prHMM->N_filtered = q.N_filtered;
      prHMM->Neff_M = (float *)calloc(prHMM->L+1, sizeof(float));
      prHMM->Neff_I = (float *)calloc(prHMM->L+1, sizeof(float));
      prHMM->Neff_D = (float *)calloc(prHMM->L+1, sizeof(float));
      prHMM->Neff_HMM = q.Neff_HMM;
      /* skip longname,name,file,fam,sfam,fold,cl */
      prHMM->lamda = q.lamda;
      prHMM->mu    = q.mu;
      
      HMMshadow rShadow = {0}; /* friend of HMM to access private members */
      rShadow.copyHMMtoShadow(q);
      
      prHMM->f     = (float **)calloc(prHMM->L+1, sizeof(float *));
      prHMM->g     = (float **)calloc(prHMM->L+1, sizeof(float *));
      prHMM->p     = (float **)calloc(prHMM->L+1, sizeof(float *));
      prHMM->tr    = (float **)calloc(prHMM->L+1, sizeof(float *));
      prHMM->linTr = (float **)calloc(prHMM->L+1, sizeof(float *));
      for (int i = 0; i < prHMM->L+1; i++){
          prHMM->f[i] = (float *)calloc(AMINOACIDS, sizeof(float));
          prHMM->g[i] = (float *)calloc(AMINOACIDS, sizeof(float));
          prHMM->p[i] = (float *)calloc(AMINOACIDS, sizeof(float));
          for (int j = 0; j < AMINOACIDS; j++){
              prHMM->f[i][j] = (float)rShadow.f[i][j];
              prHMM->g[i][j] = (float)rShadow.g[i][j];
              prHMM->p[i][j] = (float)rShadow.p[i][j];
          }
          prHMM->tr[i]    = (float *)calloc(STATE_TRANSITIONS, sizeof(float));
          prHMM->linTr[i] = (float *)calloc(STATE_TRANSITIONS, sizeof(float));
          for (int j = 0; j< STATE_TRANSITIONS; j++){
              prHMM->tr[i][j]    = (float)rShadow.tr[i][j];
              prHMM->linTr[i][j] =  fpow2(rShadow.tr[i][j]);
          }
      } /*0 <= i < prHMM->L+1 */
      /* skip trans_lin,ss_dssp,sa_dssp,ss_pred,ss_conf,Xcons */
      for (int j = 0; j < AMINOACIDS; j++){
          prHMM->pav[j] = (float)rShadow.pav[j];
      }
      /* skip pnul */
      
  } /* have read HHM from file */
  /*else if ( ((NULL != ppcProf) && (iCnt > 0) && ('\0' != ppcProf[0][0])) || 
    ( (0 != strcmp(infile,"") && (0 == iCnt) )) )*/
  else if ( (INTERN_ALN_2_HMM == iRnPtype) || (READ_ALN_2_HMM == iRnPtype) ) {
      
      if ( (INTERN_ALN_2_HMM == iRnPtype) && 
           ( (NULL == ppcProf) || (iCnt <= 0) ||  ('\0' == ppcProf[0][0]) ) ){
          printf("%s:%s:%d: was expecting internal alignment but\n"
                 "\tppcProf=%p, #seq=%d, ppcProf[0][0]=%c\n"
                 , __FUNCTION__, __FILE__, __LINE__
                 , ppcProf, iCnt, ppcProf[0][0]);
          return FAILURE;
      }
      else if ( (READ_ALN_2_HMM == iRnPtype) &&
                (0 == strcmp(infile,"")) ){
          printf("%s:%s:%d: was expecting to read alignment from file but no filename\n"
                 , __FUNCTION__, __FILE__, __LINE__);
          return FAILURE;
      }

    /*******************************/
    /*** (ii) it is an alignment ***/
    /*******************************/
    /* transfer alignment information from clustal character array
       into pali/q/t classes */
      /* 
       * NOTE that emissions in HMMer file format contain pseudo-counts.
       * HHM file format does not contain emission pseudo-counts. 
       * the structure that stores background HMM information does contain pseudo-counts
       */
    Alignment* pali;
    if (qali==NULL){ 
        pali=new Alignment(iCnt); /* FIXME: pass in iCnt to get rid of MAXSEQ */
    }
    else{
        pali=qali;
    }

    if (par.calibrate) {
      printf("\nError in %s: only HHM files can be calibrated.\n",program_name); 
      printf("Build an HHM file from your alignment with 'hhmake -i %s' and rerun hhsearch with the hhm file\n\n",infile); 
      exit(1);
    }
    
    if (v>=2 && strcmp(infile,"stdin")) cout<<infile<<" is in A2M, A3M or FASTA format\n";
    
    /* Read alignment from infile into matrix X[k][l] as ASCII 
       (and supply first line as extra argument) */ 
    //if (iCnt > 0)
    if (INTERN_ALN_2_HMM == iRnPtype){
        pali->Transfer(ppcProf, iCnt);
    }
    //else if (0 == iCnt)
    else if (READ_ALN_2_HMM == iRnPtype){
        pali->Read(inf, infile, line);
    }
    else {
        printf("%s:%s:%d: FATAL problem\n"
               "infile = (%s), #seq = %d\n"
               , __FUNCTION__, __FILE__, __LINE__
               , infile, iCnt); 
        return FAILURE;
    }
    
    /* Convert ASCII to int (0-20), throw out all insert states, 
       record their number in I[k][i] 
       and store marked sequences in name[k] and seq[k] */ 
    pali->Compress(infile);
    
    /* Sort out the nseqdis most dissimilar sequences for display 
       in the output alignments */ 
    pali->FilterForDisplay(par.max_seqid, par.coverage, par.qid,
                           par.qsc,par.nseqdis);
    
    // Filter alignment for min score per column with core query profile, defined by coverage_core and qsc_core
    //if (par.coresc>-10) pali->HomologyFilter(par.coverage_core, par.qsc_core, par.coresc);
    
    /* Remove sequences with seq. identity larger than seqid percent 
       (remove the shorter of two) */ 
    pali->N_filtered = pali->Filter(par.max_seqid, par.coverage,
                                    par.qid, par.qsc, par.Ndiff);  
    
    /* Calculate pos-specific weights, 
       AA frequencies and transitions -> f[i][a], tr[i][a] */ 
    pali->FrequenciesAndTransitions(q);
    if (v>=2 && q.Neff_HMM>11.0) 
      fprintf(stderr,"WARNING: alignment %s looks too diverse (Neff=%.1f>11). Better check it with an alignment viewer... \n",q.name,q.Neff_HMM);

    /*printf("%d %d %f %d (N,Nf,Neff,L) %s:%s:%d\n"
      , q.N_in, q.N_filtered, q.Neff_HMM, q.L, __FUNCTION__, __FILE__, __LINE__);*/

    // Add transition pseudocounts to query -> p[i][a] 
    q.AddTransitionPseudocounts();
    
    /* Generate an amino acid frequency matrix from f[i][a] 
       with full pseudocount admixture (tau=1) -> g[i][a] */ 
    q.PreparePseudocounts();
    
    /* Add amino acid pseudocounts to query:  
       p[i][a] = (1-tau)*f[i][a] + tau*g[i][a] */ 
    q.AddAminoAcidPseudocounts();      
    
    /* ****** add aligned background pseudocounts ***** */
    HMMshadow rShadowQ = {0};
    rShadowQ.copyHMMtoShadow(q);

    AddBackgroundFrequencies(rShadowQ.f, rShadowQ.p, rShadowQ.tr,
                             q.L, prHMM, 
                             q.seq, pcPrealigned, iCnt, pcRepresent);


    if (qali==NULL){
        delete(pali); pali = NULL;
    }

  } /* else if (NULL != ppcProf) // not hmmer/hhalign but alignment */

  //else if ((prHMM->L > 0) && ('\0' == ppcProf[0][0]))
  else if (INTERN_HMM_2_HMM == iRnPtype){

    /******************************************/
    /*** (iii) re-cycle old HMM information ***/
    /******************************************/

      if (prHMM->L <= 0){
          printf("%s:%s:%d: was expecting to copy HMM structure but L=%d\n"
                 , __FUNCTION__, __FILE__, __LINE__, prHMM->L);
      }

    printf("%s:%s:%d: RE-CYCLING HMM\n", __FUNCTION__, __FILE__, __LINE__);

#if 0
    q.n_display = prHMM->n_display;
    /* ignore sname*/
    for (int i = 0; i < q.n_display; i++){
      /* NOTE: In the original hhalign code the first position
         is kept open ('\0'). This makes it difficult to use
         string functions like strlen/strdup etc.
         Insert a temporary gap (.) to facilitate string operations */
      char cHead  = prHMM->seq[i][0];
      prHMM->seq[i][0] = '.';
      q.seq[i] = strdup(prHMM->seq[i]);
      q.seq[i][0] = prHMM->seq[i][0] = cHead;
    }
    q.nfirst     = prHMM->nfirst;
#else
    q.n_display  = 1; 
    q.nfirst     = 0;
    char cHead  = prHMM->seq[prHMM->nfirst][0];
    prHMM->seq[prHMM->nfirst][0] = '.';
    q.seq[0] = strdup(prHMM->seq[prHMM->nfirst]);
    q.seq[q.n_display-1][0] = prHMM->seq[prHMM->nfirst][0] = cHead;
#endif
    q.ncons      = prHMM->ncons;
    q.nss_dssp   = prHMM->nss_dssp;
    q.nsa_dssp   = prHMM->nsa_dssp;
    q.nss_pred   = prHMM->nss_pred;
    q.nss_conf   = prHMM->nss_conf;
    q.L          = prHMM->L;
    q.N_in       = prHMM->N_in;
    q.N_filtered = prHMM->N_filtered;
#define NEFF_SCORE 10 /* FIXME: Magic Number */
    /*for (int i; i < prHMM->L+1; i++){
      q.Neff_M[i] = q.Neff_I[i] = q.Neff_D[i] = NEFF_SCORE;
      }*/
    q.Neff_HMM = prHMM->Neff_HMM;
    /* skip longname,name,file,fam,sfam,fold,cl */
    q.lamda    = prHMM->lamda;
    q.mu       = prHMM->mu;

    HMMshadow rShadow = {0}; /* friend of HMM to access private members */
    rShadow.copyShadowToHMM(q, *prHMM);

  }

  if (iCnt > 0){
      if (par.forward>=1) q.Log2LinTransitionProbs(1.0);
  }

  if (NULL != inf){
      fclose(inf);
  }

  return OK;

} /*** end: ReadAndPrepare() ***/

/**
 * FreeHMMstruct()
 *
 * @brief FIXME
 *
 * @param[in,out]
 */
extern "C" void
FreeHMMstruct(hmm_light *prHMM)
{
    int i;

    if (NULL != prHMM->f){
        for (i = 0; i < prHMM->L+1; i++){
            if (NULL != prHMM->f[i]){
                free(prHMM->f[i]); prHMM->f[i] = NULL;
            }
        } /* 0 <= i < prHMM->L+1 */
        free(prHMM->f); prHMM->f = NULL;
    }
    if (NULL != prHMM->g){
        for (i = 0; i < prHMM->L+1; i++){
            if (NULL != prHMM->g[i]){
                free(prHMM->g[i]); prHMM->g[i] = NULL;
            }
        } /* 0 <= i < prHMM->L+1 */
        free(prHMM->g); prHMM->g = NULL;
    }
    if (NULL != prHMM->p){
        for (i = 0; i < prHMM->L+1; i++){
            if (NULL != prHMM->p[i]){
                free(prHMM->p[i]); prHMM->p[i] = NULL;
            }
        } /* 0 <= i < prHMM->L+1 */
        free(prHMM->p); prHMM->p = NULL;
    }
    if (NULL != prHMM->tr){
        for (i = 0; i < prHMM->L+1; i++){
            if (NULL != prHMM->tr[i]){
                free(prHMM->tr[i]); prHMM->tr[i] = NULL;
            }
        } /* 0 <= i < prHMM->L+1 */
        free(prHMM->tr); prHMM->tr = NULL;
    }
    if (NULL != prHMM->linTr){
        for (i = 0; i < prHMM->L+1; i++){
            if (NULL != prHMM->linTr[i]){
                free(prHMM->linTr[i]); prHMM->linTr[i] = NULL;
            }
        } /* 0 <= i < prHMM->L+1 */
        free(prHMM->linTr); prHMM->linTr = NULL;
    }
    
    if (NULL != prHMM->Neff_M){
        free(prHMM->Neff_M); prHMM->Neff_M = NULL;
    }
    if (NULL != prHMM->Neff_I){
        free(prHMM->Neff_I); prHMM->Neff_I = NULL;
    }
    if (NULL != prHMM->Neff_D){
        free(prHMM->Neff_D); prHMM->Neff_D = NULL;
    }

    if (NULL != prHMM->seq){
        for (i = 0; i < prHMM->n_display; i++){
            if (NULL != prHMM->seq[i]){
                free(prHMM->seq[i]); prHMM->seq[i] = NULL;
            }
        }
        free(prHMM->seq); prHMM->seq = NULL;
    }

    memset(prHMM, 0, sizeof(hmm_light));

} /*** end: FreeHMMstruct() ***/


/**
 * @brief comparisin function for ascending qsort() of doubles
 */
int 
CompDblAsc(const void *pv1, const void *pv2){

    double d1 = *(double *)pv1;
    double d2 = *(double *)pv2;

    if      (d1 > d2) { return +1; }
    else if (d1 < d2) { return -1; }
    else {              return  0; }

} /*** end: CompDblAsc() ***/


/**
 * AlnToHMM2()
 *
 * @brief convert alignment into HMM (hhmake)
 *
 * @param[out] prHMM
 * structure with pseudocounts etc
 * @param[in] pcHMM_input
 * name of file with HMM
 *
 */
extern "C" int
AlnToHMM2(hmm_light *rHMM_p, hhalign_para rHhalignPara, 
          char **ppcSeq, int iN){
    if (0 == par.M){
        SetDefaults(rHhalignPara);
        SetSubstitutionMatrix();
        par.cons = 1;
        par.M = 2;
        par.forward=2;
        par.Mgaps=100;
        const int ciGoodMeasureSeq = 10;
        int iAuxLen = strlen(ppcSeq[0]);
        par.maxResLen  = iAuxLen;
        par.maxResLen += ciGoodMeasureSeq;
        par.maxColCnt  = iAuxLen + ciGoodMeasureSeq;
        par.max_seqid=DEFAULT_FILTER;
        par.loc=0; par.mact=0;
        /* some minor parameters, affecting alignment (i think) */
        par.p = 0.0; /* minimum threshold for inclusion in hit list */
        par.Z = 100; /* minimum threshold for inclusion in hit list and alignment listing */
        par.z = 1;   /* min number of lines in hit list */
        par.B = 100; /* max number of alignments */
        par.b = 1;   /* min number of alignments */
        par.E = 1e6; /* maximum threshold for inclusion in hit list and alignment listing */
        par.altali=1;par.hitrank=0;par.showcons=1; par.showdssp=1;par.showpred=1;par.nseqdis=iN;par.cons=1;
    }

    const int ciGoodMeasure = 10;
    int iLen = strlen(ppcSeq[0]) + ciGoodMeasure;
    if (iLen > par.maxResLen){
        par.maxResLen = par.maxColCnt = iLen;
    }
    HMM rTemp(iN+2, iLen); /* temporary template */
    Alignment rTempAli(iN+2, iLen); /* temporary alignment */
    int iParCons = par.cons;

    /*printf(">>>>>>>>>>> %s:%s:%d: there are %d sequences\n", __FUNCTION__, __FILE__, __LINE__, iN);*/

    par.cons = 1;
    if (OK != ReadAndPrepare(INTERN_ALN_2_HMM, 
                             ppcSeq, iN, rHMM_p, 
                             NULL/*prealigned*/, NULL/*representative*/, NULL/*weights*/, //YES/*transfer*/,
                             (char*)("")/*in-file*/, rTemp, &rTempAli)) {
        return FAILURE;
    }
    par.cons = iParCons;

    /******/
    /*** transfer class info to struct */
    rHMM_p->n_display = rTemp.n_display;
    rHMM_p->sname = NULL;
    rHMM_p->seq = (char **)calloc((rTemp.n_display+1), sizeof(char *));

    for (int i = 0; i < rTemp.n_display; i++){
        /* NOTE: In the original hhalign code the first position
         is kept open ('\0'). This makes it difficult to use 
         string functions like strlen/strdup etc. 
         Insert a temporary gap (.) to facilitate string operations */
        char cHead  = rTemp.seq[i][0];
        rTemp.seq[i][0] = '.';
        rHMM_p->seq[i] = strdup(rTemp.seq[i]);
        rTemp.seq[i][0] = rHMM_p->seq[i][0] = cHead;
    }
    rHMM_p->ncons      = rTemp.ncons;
    rHMM_p->nfirst     = rTemp.nfirst;
    if (-1 == rHMM_p->ncons){
        /* ncons needed later for alignment of 
           representative sequence and copy of profile.
           ncons not always set (-1 default), 
           this will cause segmentation fault. 
           nfirst (probably) right index - 
           no problem if not */
        rHMM_p->ncons = rTemp.nfirst;
    }
    rHMM_p->nss_dssp   = rTemp.nss_dssp;
    rHMM_p->nsa_dssp   = rTemp.nsa_dssp;
    rHMM_p->nss_pred   = rTemp.nss_pred;
    rHMM_p->nss_conf   = rTemp.nss_conf;
    rHMM_p->L          = rTemp.L;
    rHMM_p->N_in       = rTemp.N_in;
    rHMM_p->N_filtered = rTemp.N_filtered;
#define NEFF_SCORE 10 /* FIXME: Magic Number */  //// get rid of that, FS, 2010-12-22
    rHMM_p->Neff_M = (float *)calloc(rHMM_p->L+1, sizeof(float));
    rHMM_p->Neff_I = (float *)calloc(rHMM_p->L+1, sizeof(float));
    rHMM_p->Neff_D = (float *)calloc(rHMM_p->L+1, sizeof(float));
    for (int i = 0; i < rHMM_p->L+1; i++){
        rHMM_p->Neff_M[i] = rHMM_p->Neff_I[i] = rHMM_p->Neff_D[i] = NEFF_SCORE; //// get rid of that, FS, 2010-12-22
    }
    rHMM_p->Neff_HMM = rTemp.Neff_HMM;
    /* skip longname,name,file,fam,sfam,fold,cl */
    rHMM_p->lamda = rTemp.lamda;
    rHMM_p->mu    = rTemp.mu;

    HMMshadow rShadow = {0}; /* friend of HMM to access private members */
    rShadow.copyHMMtoShadow(rTemp);

    rHMM_p->Neff_M = (float *)calloc(rHMM_p->L+1, sizeof(float));
    rHMM_p->Neff_I = (float *)calloc(rHMM_p->L+1, sizeof(float));
    rHMM_p->Neff_D = (float *)calloc(rHMM_p->L+1, sizeof(float));
    rHMM_p->f     = (float **)calloc(rHMM_p->L+1, sizeof(float *));
    rHMM_p->g     = (float **)calloc(rHMM_p->L+1, sizeof(float *));
    rHMM_p->p     = (float **)calloc(rHMM_p->L+1, sizeof(float *));
    rHMM_p->tr    = (float **)calloc(rHMM_p->L+1, sizeof(float *));
    rHMM_p->linTr = (float **)calloc(rHMM_p->L+1, sizeof(float *));

    for (int i = 0; i < rHMM_p->L+1; i++){
        rHMM_p->Neff_M[i] = (float)rShadow.Neff_M[i];
        rHMM_p->Neff_I[i] = (float)rShadow.Neff_I[i];
        rHMM_p->Neff_D[i] = (float)rShadow.Neff_D[i];
      rHMM_p->f[i] = (float *)calloc(AMINOACIDS, sizeof(float));
      rHMM_p->g[i] = (float *)calloc(AMINOACIDS, sizeof(float));
      rHMM_p->p[i] = (float *)calloc(AMINOACIDS, sizeof(float));
      for (int j = 0; j < AMINOACIDS; j++){
          rHMM_p->f[i][j] = (float)rShadow.f[i][j];
          rHMM_p->g[i][j] = (float)rShadow.g[i][j];
          rHMM_p->p[i][j] = (float)rShadow.p[i][j];
      }
      rHMM_p->tr[i]    = (float *)calloc(STATE_TRANSITIONS, sizeof(float));
      rHMM_p->linTr[i] = (float *)calloc(STATE_TRANSITIONS, sizeof(float));
      for (int j = 0; j< STATE_TRANSITIONS; j++){
          rHMM_p->tr[i][j]    = (float)rShadow.tr[i][j];
          rHMM_p->linTr[i][j] =  fpow2(rShadow.tr[i][j]);
      }
    } /*0 <= i < prHMM->L+1 */
    /* skip trans_lin,ss_dssp,sa_dssp,ss_pred,ss_conf,Xcons */
    for (int j = 0; j < AMINOACIDS; j++){
      rHMM_p->pav[j] = (float)rShadow.pav[j];
    }
    /* skip pnul */


    /******/


    rTemp.ClobberGlobal();
    rTempAli.ClobberGlobal();

    return OK;

} /*** end of AlnToHMM2() ***/


/**
 * HHMake_Wrapper()
 *
 * @brief turn alignment (from file) into HHM (HMM) on file
 *
 * @param[out] hmm_out
 * name of file with HMM info corresponding to tmp_aln
 * @param[in] tmp_aln
 * name of file with alignment, to be turned into HMM (HHM)
 *
 */
extern "C" int 
HHMake_Wrapper(char *tmp_aln, char *hmm_out)
{

    HMM rTemp; /* temporary template */
    Alignment rTempAli; /* temporary alignment */
    hmm_light *rHMM_p = NULL;

    /** Note:
        this is a wrapper for a stand-alone program hhmake.
        hhmake uses a different set of parameters than hhalign.
        However, all parameters are GLOBAL. 
        So at this point we register the hhalign settings, 
        reset them to hhmake settings and set them back 
        at the end of the function
     */

    /* save settings from hhalign */
    int iParShowcons=par.showcons;
    int iParAppend=par.append;
    int iParNseqdis=par.nseqdis;
    int iParMark=par.mark;
    int iParMaxSeqid=par.max_seqid;
    int iParQid=par.qid;
    double dParQsc=par.qsc;
    int iParCoverage=par.coverage;
    int iParNdiff=par.Ndiff;
    int iParCoverageCore=par.coverage_core;
    double dParQscCore=par.qsc_core;
    double dParCoresc=par.coresc;
    int iParM=par.M;
    int iParMgaps=par.Mgaps;
    int iParMatrix=par.matrix;
    int iParPcm=par.pcm;
    double dParPcw=par.pcw;
    double dParGapb=par.gapb;
    int iParWg=par.wg;

    /* these are settings suitable for hhmake */
    par.showcons=1;              // write consensus sequence into hhm file
    par.append=0;                // overwrite output file                
    par.nseqdis=10;              // maximum number of query or template sequences to be recoreded in HMM and diplayed in output alignments 
    par.mark=0;                  // 1: only marked sequences (or first) get displayed; 0: most divergent ones get displayed 
    par.max_seqid=90;            // default for maximum sequence identity threshold
    par.qid=0;                   // default for maximum sequence identity threshold
    par.qsc=-20.0f;              // default for minimum score per column with query
    par.coverage=0;              // default for minimum coverage threshold         
    par.Ndiff=100;               // pick Ndiff most different sequences from alignment
    par.coverage_core=30;        // Minimum coverage for sequences in core alignment  
    par.qsc_core=0.3f;           // Minimum score per column with core alignment (HMM)
    par.coresc=-20.0f;           // Minimum score per column with core alignment (HMM)
    par.M=1;                     // match state assignment is by A2M/A3M              
    par.Mgaps=50;                // above this percentage of gaps, columns are assigned to insert states
    par.matrix=0;                // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50 3: BLOSUM62            
    par.pcm=0;                   // no amino acid and transition pseudocounts added                     
    par.pcw=0;                   // wc>0 weighs columns according to their intra-clomun similarity      
    par.gapb=0.0;                // default values for transition pseudocounts; 0.0: add no transition pseudocounts!
    par.wg=0;                    // 0: use local sequence weights   1: use local ones                               

    
    if (OK != ReadAndPrepare(READ_ALN_2_HMM, 
                             NULL, 0, rHMM_p, NULL, NULL, NULL,
                             tmp_aln, rTemp, &rTempAli)) {
        return FAILURE;
    }

    rTemp.WriteToFile(hmm_out);



    /* restore settings from hhalign */
    par.showcons=iParShowcons;
    par.append=iParAppend;
    par.nseqdis=iParNseqdis;
    par.mark=iParMark;
    par.max_seqid=iParMaxSeqid;
    par.qid=iParQid;
    par.qsc=dParQsc;
    par.coverage=iParCoverage;
    par.Ndiff=iParNdiff;
    par.coverage_core=iParCoverageCore;
    par.qsc_core=dParQscCore;
    par.coresc=dParCoresc;
    par.M=iParM;
    par.Mgaps=iParMgaps;
    par.matrix=iParMatrix;
    par.pcm=iParPcm;
    par.pcw=dParPcw;
    par.gapb=dParGapb;
    par.wg=iParWg;


    /* prepare exit */
    rTemp.ClobberGlobal();
    rTempAli.ClobberGlobal();

    return OK;
}

/**
 * readHMMWrapper()
 *
 * @brief read HMM from file, copy (relevant) info into struct
 *
 * @param[out] prHMM
 * structure with pseudocounts etc, probably uninitialised on entry
 * @param[in] pcHMM_input
 * name of file with HMM
 *
 */
extern "C" int
readHMMWrapper(hmm_light *rHMM_p, 
                 char *pcHMM_input)
{
    par.maxResLen = 15002;
    HMM rTemp(1000,par.maxResLen); /* temporary template */
    Alignment rTempAli; /* temporary alignment */
    /* AW changed init from {0} to 0 because it failed to compile with
     * gcc 4.3.3 with the following error:
     * error: braces around initializer for non-aggregate type 
     */
    /* FS taken out initialiser alltogether */    

    /* 0th arg of RnP is the type of RnP, ie, 
       enum {INTERN_ALN_2_HMM = 0, READ_ALN_2_HMM, READ_HMM_2_HMM};*/
    /* 1st arg of ReadAndPrepare() is ppcSeqs, 2nd is #seq */
    /* FIXME: at this stage the 3rd arg, rHMM_p, only acts as a dummy.
       this is rather silly, as it is the actual struct that will 
       carry the HMM info at the end. 
       If we write it already in ReadAndPrepare() then we could 
       dispense with all this friend-class nonsense */
    if (OK != ReadAndPrepare(READ_HMM_2_HMM, 
                             NULL, 0, rHMM_p, NULL, NULL, NULL, 
                             pcHMM_input, rTemp, &rTempAli)) {
        return FAILURE;
    }

    /* an important piece of information I want to get out of here 
       is the consenssus sequence. there are however certain 
       Pfam HMMs that don't trigger consensus calculation.
       at the moment I (FS) don't understand why this is 
       (or rather why this is not). The proper place to do this 
       should be inside ReadAndPrepare():ReadHMMer(), but 
       I have not quite penetrated the logic there. 
       Therefore I try to catch this condition at this point (here) 
       and rectify it. 
     */
    if (-1 == rHMM_p->ncons){
        int i, iA;
        rHMM_p->ncons = rHMM_p->nfirst;

        for (i = 0; i < rHMM_p->L; i++){
            double dPmax = 0.00;
            int iAmax = -1;
            for (iA = 0; iA < AMINOACIDS; iA++){
                if (rHMM_p->f[i][iA] > dPmax){
                    iAmax = iA;
                    dPmax = rHMM_p->f[i][iA];
                }
            } /* (0 <= iA < AMINOACIDS) */
            rHMM_p->seq[rHMM_p->ncons][i] = i2aa(iAmax);
        } /* (0 <= i < rHMM_p->L) */
        rHMM_p->seq[rHMM_p->ncons][i] = '\0'; /* FS, r291 -> */
    } /* ncons not set */


    rTemp.ClobberGlobal();
    rTempAli.ClobberGlobal();

  return OK;

}  /*** end: readHMMWrapper() ***/






/////////////////////////////////////////////////////////////////////////////
/**
 * @brief Do precalculations for q and t to prepare comparison
 */
void 
PrepareTemplate(HMM& q, HMM& t, int format)
{
  if (format==0) // HHM format
    {
      // Add transition pseudocounts to template
      t.AddTransitionPseudocounts();

      // Generate an amino acid frequency matrix from t.f[i][a] with max pseudocounts (tau=1) -> t.g[i][a] 
      t.PreparePseudocounts();

      // Add amino acid pseudocounts to template: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
      t.AddAminoAcidPseudocounts();
    }
  else // HHMER format
    {
      // Don't add transition pseudocounts to template
      // t.AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg, par.gaph, par.gapi, 0.0);
      
      // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a] 
      t.PreparePseudocounts();

      // DON'T ADD amino acid pseudocounts to temlate: pcm=0!  t.p[i][a] = t.f[i][a]
      t.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
     }

  // Modify transition probabilities to include SS-dependent penalties
  if (par.ssgap) t.UseSecStrucDependentGapPenalties();

  if (par.forward>=1) t.Log2LinTransitionProbs(1.0);  

  // Factor Null model into HMM t
  // ATTENTION! t.p[i][a] is divided by pnul[a] (for reasons of efficiency) => do not reuse t.p
  t.IncludeNullModelInHMM(q,t);  // Can go BEFORE the loop if not dependent on template
  
  return;
}


/***    end of hhfunc-C.h   ***/
