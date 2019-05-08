/* hpm_scoreops.c */

#include "easel.h"
#include "hmmer.h"
#include "p7_config.h"
#include "esl_sq.h"
#include "hpm.h"
#include "hpm_scoreops.h"

int hpm_scoreops_CalculateHamiltonian(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_hsc, float *ret_esc) {
   int   z;         /* index for trace elements                   */
   int   y;         /* index for trace elements                   */
   int   i;         /* match state index                          */
   int   j;         /* match state index                          */
   int   a;         /* residue index                              */
   int   b;         /* residue index                              */
   int   idx;       /* index for potts parameters                 */
   float hsc = 0;   /* Hamiltonian contribution from h_i's        */
   float esc = 0;   /* Hamiltonian contribution from e_ij's       */


   int K = hpm->abc->K;

   for (z = 0; z < tr->N; z++) {

      //fprintf(stdout, "%d\n", z);
      /* check if we are in a match or delete state */
      if (tr->st[z] == 2 || tr->st[z] == 6) {
         i = tr->k[z];

         /* we have a match position */
         if (tr->st[z] == 2) {
            a = dsq[tr->i[z]];
            if (a > K) a = K;

            //fprintf(stdout, "%d\n", a);
         }

         /* we have a delete position */
         else if (tr->st[z] == 6) {
            a = K;
         }
         hsc += hpm->h[i][a];

         /* now add e_ij terms to pseudo-energy */
         for (y = z+1; y < tr->N; y++) {

            /* check if we are in a match position */
            if (tr->st[y] == 2 || tr->st[y] == 6) {
               j = tr->k[y];

               /* we have a match state */
               if (tr->st[y] == 2) {
                  b = dsq[tr->i[y]];
                  if (b > K) b = K;
               }

               /* we have a delete state */
               else if (tr->st[y] == 6) {
                  b = K;
               }

               idx = IDX(a,b,K+1);
               //fprintf(stdout, "i=%d, j=%d, a=%d, b=%d, idx=%d\n", i, j, a, b, idx);
               esc += hpm->e[i][j][idx];
            }
         }
      }
   }

   *ret_hsc = hsc;
   *ret_esc = esc;
   return eslOK;
}


int hpm_scoreops_CalculateHamiltonianSingleSite(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, int i, float *ret_Ei) {
   int   y;           /* index for trace elements                   */
   int   z;           /* index for trace elements                   */
   int   j;           /* match state index                          */
   int   a;           /* residue index                              */
   int   b;           /* residue index                              */
   int   idx;         /* index for potts parameters                 */
   float Ei = 0.0;   /* Hamiltonian contribution for site i        */
   int K = hpm->abc->K;

	//fprintf(stdout, "\t in hpm_scoreops_CalculateHamiltonianSingleSite()\n");

	/* figure out which node and residue we are dealing with */
	for (y = 0; y < tr->N; y++) {

		//fprintf(stdout, "\t\ty: %d\n", y);
		if (tr->st[y] == p7T_MG || tr->st[y] == p7T_MG) {
			if (tr->k[y] == i) {
				z = y;
				a = dsq[tr->i[z]];
				//fprintf(stdout, "\t\t\tWe've found our node! z = %d k = %d a =%d\n",
				//		  z, tr->k[y], a);
				break;
			}
		}
	}

	/* get h_i term */
	Ei += hpm->h[i][a];
	for (y = 0; y < tr->N; y++) {
		if (z == y) continue;
		//fprintf(stdout, "\t\ty: %d\n", y);

		if (tr->st[y] == p7T_MG || tr->st[y] == p7T_DG) {

			j = tr->k[y];

			/* we have a match state */
			if (tr->st[y] == p7T_MG) {
				b = dsq[tr->i[y]];
				if (b > K) b = K;
			}

			/* we have a delete state */
			else if (tr->st[y] == p7T_DG) {
				b = K;
			}

			idx = IDX(a,b,K+1);
			//fprintf(stdout, "\t\t\ti=%d, j=%d, a=%d, b=%d, idx=%d, eij[a][b]=%.4f\n",
			//		  i, j, a, b, idx,hpm->e[i][j][idx]);
         Ei += hpm->e[i][j][idx];
		}
	}

	/* loop over trace to get e_ij terms */

   *ret_Ei = Ei;

   return eslOK;
}




int hpm_scoreops_ScoreInsertEmissions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_isc)
{
	int     z;              /* index for trace elements   */
	int     i;              /* match state index          */
	int     a;              /* residue index            	*/
	float   isc = 0.0;

	for (z = 0; z < tr->N; z++) {

		/* we have an insert position */
		/* IG, N, or C states */
		if (tr->st[z] == p7T_IG || (tr->st[z] == p7T_N && tr->i[z] > 0 ) || ((tr->st[z] == p7T_C && tr->i[z] > 0 )) ) {

			/* get node index */
			i = tr->k[z];
			/* get emitted residue */
			a = dsq[tr->i[z]];
			/* update insert emission log prob */
			isc += log(hpm->ins[i][a]);
		}
	}

	*ret_isc = isc;
	return eslOK;
}

int hpm_scoreops_ScoreNullEmissions(HPM *hpm, ESL_SQ *sq, float *ret_nesc) {
   int    i;                 /* match state index             */
   int    a;                 /* residue index                 */
   float  nesc = 0.0;        /* log prob of insert emissions  */
   int    K = hpm->abc->K;   /* alphabet size                 */

   for (i = 1; i < sq->n+1; i++) {
      a = sq->dsq[i];
      /* treat all degenerate residues as first letter in abc  */
      /* insert emissions cancel in log odds score anyways     */
      /* this is possibly a half baked idea, 11/13/2018        */
      if (a > K) a = 0;
         nesc += log(hpm->ins[0][a]);

   }

   *ret_nesc = nesc;
   return eslOK;
}


int hpm_scoreops_ScoreNullMatchEmissions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_nmesc) {
   int     z;            /* index for trace elements       */
   int     i;            /* match state index                 */
   int     a;            /* residue index                  */
   float   nmesc = 0.0;  /* log prob of insert emissions     */

   int K = hpm->abc->K;

   for (z = 0; z < tr->N; z++) {
      /* we have a match state */
      if (tr->st[z] == 2) {
         i = tr->k[z];
         a = dsq[tr->i[z]];
         /* treat all degenerate residues as first letter in abc  */
         /* insert emissions cancel in log odds score anyways     */
         /* this is possibly a half baked idea, 11/13/2018        */
         if (a > K) a = 0;
         nmesc += log(hpm->ins[i][a]);
      }
   }

   *ret_nmesc = nmesc;

   return eslOK;
}


int hpm_scoreops_ScoreTransitions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, int L, float *ret_tsc) {
   int   z;              /* index for trace elements               */
   int   st;
   int   stprev  = -1;   /* state id index                         */
   int   i;              /* match state index                      */
   int   iprev;          /* match state index                      */
   float tsc     = 0.0;  /* log of transition probs                */
   float tsc_NN  = 0.0;  /* log of N->N and C->C transition prob   */
   float tsc_NB  = 0.0;  /* log of N->B and C->T transition prob   */

   /* calculate N->N and C->C transition probability, in log space */
   tsc_NN = log(L) - log(L+2);
   /* calculate N->B and C->T transition probability, in log space */
   tsc_NB = log(2) - log(L+2);

   /* loop over trace positions for this seq */
   for (z = 0; z < tr->N; z++) {
      /* state type */
      st = tr->st[z];
      /* node in model */
      i = tr->k[z];

      /* handle transitions into N state */
      if (st == 8) {
         /*ignore S->N transitions, they have prob 1 */

         /* handle N->N transitions */
         if (stprev == 8)                     tsc = tsc + tsc_NN;
      }

      /* handle transitions into B state */
      else if (st == 9) {
         /* handle N->B transitons */
         if (stprev == 8)                     tsc += tsc_NB;
         /* no other transitions into B state possible */
      }

      /* handle transitions into hybrid match-del states */
      else if (st == 2 || st == 6) {
         /* handle M->M transitions */
         if      (stprev == 2 || stprev == 6) tsc += log(hpm->t[iprev][HPM_MM]);
         /* handle I->M transitions */
         else if (stprev == 4)                tsc += log(hpm->t[iprev][HPM_IM]);
         /* ignore B->M transitions, they have prob 1 */
      }
      /* handle transitions into I states */
      else if (st == 4) {
         /* handle M->I transitions */
         if      (stprev == 2 || stprev == 6) tsc += log(hpm->t[iprev][HPM_MI]);
         /* handle I->I transitions */
         else if (stprev == 4)                tsc += log(hpm->t[iprev][HPM_II]);
      }

      /* handle transitions into C state */
      else if (st == 13) {
         if (stprev == 13)                    tsc += tsc_NN;
         /* ignore E->C transitions, they have prob 1 */
      }

      /*handle transitions into T state */
      else if (st == 15) {
         if (stprev == 13)                    tsc += tsc_NB;
         /* no other possible transitions to T state */
      }

      stprev = st;
      iprev = i;
   }

   *ret_tsc = tsc;

   return eslOK;
}

int hpm_scoreops_ScoreNullTransitions(P7_BG  *bg, ESL_SQ *sq, float *ret_ntsc) {
   float ntsc = 0.0;  /* log of null transition probs */

   p7_bg_SetLength     (bg, sq->n);
   p7_bg_NullOne(bg, sq->dsq, sq->n, &ntsc);

   *ret_ntsc = ntsc;

   return eslOK;

}
