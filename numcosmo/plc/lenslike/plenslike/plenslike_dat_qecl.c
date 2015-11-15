#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "plenslike.h"

int load_plenslike_dat_qecl( plenslike_dat_qecl *dat, char *tfname ) {
  char line[32767];
  FILE *tf;
  int i, j, l, err, tmp;
  double dmp;
  char cmp1, cmp2;
  
  tf = fopen(tfname, "r");
  if (tf == NULL) {
    perr("Error opening data file.");
    return -1;
  }

  // bypass comments
  while (fgets(line, sizeof line, tf)) {
    if (*line == '#') {
      continue;
    } else {
      break;
    }
  }

  // read size header
  err = sscanf(line, "%d", &dat->nbins);
  err = fscanf(tf,   "%d %d", &dat->lmaxcmb, &dat->lmaxphi);
  err = fscanf(tf,   "%d %d", &dat->nqe, &dat->nx);

  // allocate memory for bins and spectra.
  dat->bin_lmins     = malloc( dat->nbins * sizeof(int) );
  dat->bin_lmaxs     = malloc( dat->nbins * sizeof(int) );
  dat->bin_vals      = malloc( dat->nbins * sizeof(double) );
  dat->mat_sigma     = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->mat_sigma_inv = malloc( dat->nbins * dat->nbins * sizeof(double) );

  dat->cltt_fid      = malloc( (dat->lmaxcmb+1) * sizeof(double) );
  dat->clee_fid      = malloc( (dat->lmaxcmb+1) * sizeof(double) );
  dat->clte_fid      = malloc( (dat->lmaxcmb+1) * sizeof(double) );

  dat->clpp_fid      = malloc( (dat->lmaxphi+1) * sizeof(double) );
  dat->qlpp_fid      = malloc( (dat->lmaxphi+1) * sizeof(double) );  
  dat->vlpp_inv      = malloc( (dat->lmaxphi+1) * sizeof(double) );
  dat->qlpp_inv      = malloc( (dat->lmaxphi+1) * sizeof(double) );

  dat->qes           = malloc( dat->nqe * sizeof(qest *) );
  dat->qefs          = malloc( dat->nqe * sizeof(int) );
  dat->qe12          = malloc( dat->nx  * sizeof(int) );
  dat->qe34          = malloc( dat->nx  * sizeof(int) );
  
  // read bin info
  for (i=0; i<dat->nbins; i++) {
    err = fscanf(tf, "%d %d %d %lf", &tmp, &dat->bin_lmins[i], &dat->bin_lmaxs[i], &dat->bin_vals[i]);

    if ( ( err != 4 ) || ( tmp != i ) ) {
	perr("Error reading bin info.");
	free_plenslike_dat_qecl(dat);
	return -1;
    }
  }

  // read sigma matrix
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      err = fscanf(tf, "%lf", &dat->mat_sigma[i*dat->nbins + j]);
      if (err != 1) {
	perr("Error reading sigma matrix.");
	free_plenslike_dat_qecl(dat);
	return -1;
      }
    }
  }

  // read sigma inv matrix
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      err = fscanf(tf, "%lf", &dat->mat_sigma_inv[i*dat->nbins + j]);
      if (err != 1) {
	perr("Error reading sigma^{-1} matrix.");
	free_plenslike_dat_qecl(dat);
	return -1;
      }
    }
  }

  // read cmb spectra
  for (l=0; l<=dat->lmaxcmb; l++) {
    err = fscanf(tf, "%lf %lf %lf %lf", &dmp, &dat->cltt_fid[l], &dat->clee_fid[l], &dat->clte_fid[l]);
    if ( (err != 4) || ( dmp != l) ) {
	perr("Error reading CMB spectra.");
	free_plenslike_dat_qecl(dat);
	return -1;
    }
  }

  // read phi spectra
  for (l=0; l<=dat->lmaxphi; l++) {
    err = fscanf(tf, "%lf %lf %lf %lf %lf", &dmp, &dat->clpp_fid[l], &dat->qlpp_fid[l], &dat->vlpp_inv[l], &dat->qlpp_inv[l]);
     if ( (err != 5) || ( dmp != l) ) {
	perr("Error reading PHI spectra.");
	free_plenslike_dat_qecl(dat);
	return -1;
    }
  }

  // read quadratic estimators
  for (i=0; i < dat->nqe; i++) {
    dat->qes[i] = malloc( sizeof(qest) );
    dat->qes[i]->s12L = NULL;
    dat->qes[i]->w12L = NULL;
    
    err = fscanf(tf, "%d", &dat->qes[i]->ntrm);
    err = fscanf(tf, "%d", &dat->qes[i]->lmax);

    dat->qes[i]->s12L = malloc( dat->qes[i]->ntrm*sizeof(int *) );
    for (j=0; j < dat->qes[i]->ntrm; j++) { dat->qes[i]->s12L[j] = NULL; };

    for (j=0; j < dat->qes[i]->ntrm; j++) {
      dat->qes[i]->s12L[j] = malloc( 3*sizeof(int) );
      err = fscanf(tf, "%d %d %d", &dat->qes[i]->s12L[j][0], &dat->qes[i]->s12L[j][1], &dat->qes[i]->s12L[j][2]);
      if (err != 3) {
	 perr("Error reading spins for quadratic estimator %d.", i);
	 free_plenslike_dat_qecl(dat);
	 return -1;
      }
    }
  
    dat->qes[i]->w12L = malloc( dat->qes[i]->ntrm*sizeof(double **) );
    for (j=0; j < dat->qes[i]->ntrm; j++) { dat->qes[i]->w12L[j] = NULL; };

    for (j=0; j < dat->qes[i]->ntrm; j++) {
      dat->qes[i]->w12L[j]    = malloc( 3*sizeof(double *) );
      dat->qes[i]->w12L[j][0] = malloc( (dat->qes[i]->lmax+1)*sizeof(double) );
      dat->qes[i]->w12L[j][1] = malloc( (dat->qes[i]->lmax+1)*sizeof(double) );
      dat->qes[i]->w12L[j][2] = malloc( (dat->qes[i]->lmax+1)*sizeof(double) );
      
      for (l=0; l <= dat->qes[i]->lmax; l++) {
	err = fscanf(tf, "%lf %lf %lf %lf", &dmp, &dat->qes[i]->w12L[j][0][l], &dat->qes[i]->w12L[j][1][l], &dat->qes[i]->w12L[j][2][l] );
	assert( err == 4 );
	assert( dmp == l );
	if ( (err != 4) || (dmp != l) ) {
	  perr("Error reading weights for quadratic estimator %d. Term = %d,  l = %d.", i, j, l);
	  free_plenslike_dat_qecl(dat);
	  return -1;
	}
      }
    }
  }
  
  // read fields
  for (i=0; i < dat->nqe; i++) {
    err = fscanf(tf, "%d %c %c", &tmp, &cmp1, &cmp2);
    if ( (tmp != i) || (err != 3) ) {
      perr("Error reading fields for quadratic estimator %d.", i);
      free_plenslike_dat_qecl(dat);
      return -1;
    }
    
    switch(cmp1) 
      {
      case 'T':
	dat->qefs[i] = 0*3;
	break;
      case 'E':
	dat->qefs[i] = 1*3;
	break;
      case 'B':
	dat->qefs[i] = 2*3;
	break;
      default:
	perr("Unphysical field for quadratic estimator %d.", i);
	free_plenslike_dat_qecl(dat);
	return -1;
      }

    switch (cmp2) 
      {
      case 'T':
	dat->qefs[i] += 0;
	break;
      case 'E':
	dat->qefs[i] += 1;
	break;
      case 'B':
	dat->qefs[i] += 2;
	break;
      default:
	perr("Unphysical field for quadratic estimator %d.", i);
	free_plenslike_dat_qecl(dat);
	return -1;
      }
  }

  // read qe12, qe34
  for (i=0; i < dat->nx; i++) {
    err = fscanf(tf, "%d %d %d", &tmp, &dat->qe12[i], &dat->qe34[i]);
    if ( (tmp != i) || (err != 3) ) {
      perr("Error reading cross-indices for term %d.", i);
      free_plenslike_dat_qecl(dat);
      return -1;
    }
  }

  return 0;
};

void free_plenslike_dat_qecl( plenslike_dat_qecl *dat ) {
  int i;

  free(dat->bin_lmins);
  free(dat->bin_lmaxs);
  free(dat->bin_vals);
  free(dat->mat_sigma);
  free(dat->mat_sigma_inv);
  free(dat->cltt_fid);
  free(dat->clee_fid);
  free(dat->clte_fid);

  free(dat->clpp_fid);
  free(dat->qlpp_fid);
  free(dat->vlpp_inv);
  free(dat->qlpp_inv);

  for (i=0; i < dat->nqe; i++) {
    if (dat->qes[i] != NULL) {
      free_qe(dat->qes[i]);
      free(dat->qes[i]);
    }
  }
  free(dat->qes);
  free(dat->qefs);
  free(dat->qe12);
  free(dat->qe34);
};

void fill_plenslike_qecl_bins( plenslike_dat_qecl *dat, double *bins, double *clpp ) {
  int i, l;
  double num, den;

  for (i=0; i<dat->nbins; i++) {
    num = 0; den=0;
    for (l=dat->bin_lmins[i]; l<=dat->bin_lmaxs[i]; l++) {
      num += clpp[l] * dat->clpp_fid[l] * dat->vlpp_inv[l];
      den += dat->clpp_fid[l] * dat->clpp_fid[l] * dat->vlpp_inv[l];
    }
    bins[i] = num/den;
  }
}

void fill_qecl_resp_pp_cls( int lmax, double *rl, plenslike_dat_qecl *dat, double *cltt, double *clee, double *clte ) {
  qest  **qls = malloc( 9*sizeof(qest *) );

  init_qls( qls, dat->lmaxcmb, cltt, clee, clte );

  fill_qecl_resp_pp_qls( lmax, rl, dat, qls );

  free_qls(qls);
  free(qls);
}

void fill_qecl_resp_pp_qls( int lmax, double *rl, plenslike_dat_qecl *dat, qest **qls ) {
  int i, l;
  double *fl, *respqe;

  fl = malloc( (dat->lmaxcmb+1) * sizeof(double) ); for (l=0; l<=dat->lmaxcmb; l++) { fl[l] = 1.0; }
  
  respqe = malloc( (lmax+1) * dat->nqe * sizeof(double) ); 
  memset( respqe, 0, (lmax+1) * dat->nqe * sizeof(double) );

  for (i=0; i<dat->nqe; i++) {
    fill_qe_resp(lmax, &respqe[(lmax+1)*i ], 
		 dat->qes[i],
		 qls[ dat->qefs[i] ],
		 fl, dat->qes[i]->lmax, 
		 fl, dat->qes[i]->lmax);
  }

  memset( rl, 0, (lmax+1) * sizeof(double) );
  for (i=0; i<dat->nx; i++) {
    for (l=0; l<=lmax; l++) {
      rl[l] += respqe[(lmax+1)*dat->qe12[i] + l] * respqe[(lmax+1)*dat->qe34[i] + l];
    }
  }

  free(respqe);
  free(fl);
};

double calc_plenslike_qecl( plenslike_dat_qecl *dat, double *clpp ) {
  int i, j;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );

  fill_plenslike_qecl_bins(dat, bins, clpp);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  free(bins);

  return -0.5*ret;
};

double calc_plenslike_qecl_renorm( plenslike_dat_qecl *dat, double *clpp, double *cltt, double *clee, double *clte ) {
  int i, j, l;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );
  double *clpp_renorm = malloc( (dat->lmaxphi+1) * sizeof(double) );
  double *qc_resp = malloc( (dat->lmaxphi+1) * sizeof(double) );
  
  fill_qecl_resp_pp_cls( dat->lmaxphi, qc_resp, dat, cltt, clee, clte );

  for (l=0; l<=dat->lmaxphi; l++) {
    clpp_renorm[l] = qc_resp[l] * dat->qlpp_inv[l] * clpp[l];
  }

  fill_plenslike_qecl_bins(dat, bins, clpp_renorm);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  free(bins);
  free(clpp_renorm);
  free(qc_resp);

  return -0.5*ret;
};
