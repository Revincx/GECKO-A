/* Database module for GECKO-A */
#ifndef DATABASE_H
#define DATABASE_H

#include "keyparameter.h"
#include "references.h"
#include <stdbool.h>

/* Rate constants for VOC+oxidant */
#define MXKDB 1000  /* max # of rate constant in a given database */

/* VOC+OH database */
extern int nkohdb;                              /* # of species in the database */
extern char kohdb_chem[MXKDB][MXLFO+1];         /* formula of the species */
extern double kohdb_298[MXKDB];                 /* rate constant @ 298 K */
extern double kohdb_arr[MXKDB][3];              /* arrhenius parameter for the rate constant */
extern char kohdb_com[MXKDB][3][MXLCOD+1];      /* comment's code for the rate */

/* VOC+O3 database */
extern int nko3db;                              /* # of species in the database */
extern char ko3db_chem[MXKDB][MXLFO+1];         /* formula of the species */
extern double ko3db_298[MXKDB];                 /* rate constant @ 298 K */
extern double ko3db_arr[MXKDB][3];              /* arrhenius parameter for the rate constant */
extern char ko3db_com[MXKDB][3][MXLCOD+1];      /* comment's code for the rate */

/* VOC+NO3 database */
extern int nkno3db;                             /* # of species in the database */
extern char kno3db_chem[MXKDB][MXLFO+1];        /* formula of the species */
extern double kno3db_298[MXKDB];                /* rate constant @ 298 K */
extern double kno3db_arr[MXKDB][3];             /* arrhenius parameter for the rate constant */
extern char kno3db_com[MXKDB][3][MXLCOD+1];     /* comment's code for the rate */

/* Known mechanism for VOC+oxidant ("kw" data base) */
#define MXKWR 110  /* max # of known reaction in a given database */

/* VOC+OH mechanism */
extern int nkwoh;                               /* # of species in the database */
extern char kwoh_rct[MXKWR][MXLFO+1];           /* formula of the VOC reacting with OH */
extern int nkwoh_pd[MXKWR];                     /* # of channel involved in mechanism */
extern double kwoh_yld[MXKWR][MXNR];            /* yield of each channel */
extern char kwoh_pd[MXKWR][MXNR][MXLFO+1];      /* formula of the main product per channel */
extern char kwoh_copd[MXKWR][MXNR][MXCOPD][MXLCO+1]; /* coproducts per channel */
extern char kwoh_com[MXKWR][3][MXLCOD+1];       /* comment's code for the rate */

/* VOC+O3 mechanism */
extern int nkwo3;                               /* # of species in the database */
extern char kwo3_rct[MXKWR][MXLFO+1];           /* formula of the VOC reacting with O3 */
extern int nkwo3_pd[MXKWR];                     /* # of channel involved in mechanism */
extern double kwo3_yld[MXKWR][MXNR];            /* yield of each channel */
extern char kwo3_pd[MXKWR][MXNR][MXLFO+1];      /* formula of the main product per channel */
extern char kwo3_copd[MXKWR][MXNR][MXCOPD][MXLCO+1]; /* coproducts per channel */
extern char kwo3_com[MXKWR][3][MXLCOD+1];       /* comment's code for the rate */

/* VOC+NO3 mechanism */
extern int nkwno3;                              /* # of species in the database */
extern char kwno3_rct[MXKWR][MXLFO+1];          /* formula of the VOC reacting with NO3 */
extern int nkwno3_pd[MXKWR];                    /* # of channel involved in mechanism */
extern double kwno3_yld[MXKWR][MXNR];           /* yield of each channel */
extern char kwno3_pd[MXKWR][MXNR][MXLFO+1];     /* formula of the main product per channel */
extern char kwno3_copd[MXKWR][MXNR][MXCOPD][MXLCO+1]; /* coproducts per channel */
extern char kwno3_com[MXKWR][3][MXLCOD+1];      /* comment's code for the rate */

/* Known mechanism for radicals (RO2, RCOO2, RO and criegee) */

/* RO2 chemistry */
extern int nkwro2;                              /* # of entry (reaction) in the database */
extern double kwro2_arrh[MXKWR][3];             /* arrhenius parameter for the rate constant */
extern double kwro2_stoi[MXKWR][4];             /* stoi. coef. of the products */
extern char kwro2_rct[MXKWR][2][MXLFO+1];       /* formula of the reactants */
extern char kwro2_prd[MXKWR][4][MXLFO+1];       /* formula of the products */
extern char kwro2_com[MXKWR][4][MXLCOD+1];      /* comment's code for the ro2 reaction */

/* RCOO2 chemistry */
extern int nkwrco3;                             /* # of entry (reaction) in the database */
extern double kwrco3_arrh[MXKWR][3];            /* arrhenius parameter for the rate constant */
extern double kwrco3_stoi[MXKWR][4];            /* stoi. coef. of the products */
extern char kwrco3_rct[MXKWR][2][MXLFO+1];      /* formula of the reactants */
extern char kwrco3_prd[MXKWR][4][MXLFO+1];      /* formula of the products */
extern char kwrco3_com[MXKWR][4][MXLCOD+1];     /* comment's code for the rcoo2 reaction */

/* RO chemistry */
extern int nkwro;                               /* # of entry (reaction) in the database */
extern double kwro_arrh[MXKWR][3];              /* arrhenius parameter for the rate constant */
extern double kwro_stoi[MXKWR][4];              /* stoi. coef. of the products */
extern char kwro_rct[MXKWR][2][MXLFO+1];        /* formula of the reactants */
extern char kwro_prd[MXKWR][4][MXLFO+1];        /* formula of the products */
extern char kwro_com[MXKWR][4][MXLCOD+1];       /* comment's code for the ro reaction */

/* Criegee chemistry */
extern int nkwcri;                              /* # of entry (reaction) in the database */
extern double kwcri_arrh[MXKWR][3];             /* arrhenius parameter for the rate constant */
extern double kwcri_stoi[MXKWR][4];             /* stoi. coef. of the products */
extern char kwcri_rct[MXKWR][2][MXLFO+1];       /* formula of the reactants */
extern char kwcri_prd[MXKWR][4][MXLFO+1];       /* formula of the products */
extern char kwcri_com[MXKWR][3][MXLCOD+1];      /* comment's code for the ro reaction */

/* photolysis labels & known photolysis reactions */
extern int njdat;                               /* # of known photo. reactions */
extern char jchem[MXKWR][MXLFO+1];              /* formula of the species being photolyzed */
extern char jprod[MXKWR][2][MXLFO+1];           /* the two main photodissociation fragments */
extern int jlabel[MXKWR];                       /* ID number of photolytic reaction i */
extern char coprodj[MXKWR][MXLCO+1];            /* additional inorganic coproduct */
extern char j_com[MXKWR][3][MXLCOD+1];          /* comment's code for the photolysis reactions */
extern int nj40;                                /* # of J40 data in database */
extern int jlab40[MXKWR];                       /* label for which J40 is provided */
extern double j40[MXKWR];                       /* J40 values for various labels */

/* Benson groups database */
#define MXBG 500    /* max # of benson group */
#define LBG 24      /* max length of a benson group (string) */
extern int nbson;                               /* # of group in the database */
extern double bsonval[MXBG];                    /* values of Benson group */
extern char bsongrp[MXBG][LBG+1];               /* Benson group string */

/* Simpol SAR for vapor pressure */
extern double p0simpgrp298;                     /* constant term for psat @ 298K */
extern double psimpgrp298[30];                  /* grp contribution for psat @ 298K */
extern double h0simpgrp298;                     /* constant term for heat of vaporisation @ 298K */
extern double hsimpgrp298[30];                  /* grp contribution for heat of vaporisation @ 298K */

/* Nannoolal SAR for vapor pressure and Tb */
extern double nanweight[219][2];                /* grp contribution. 1st column is Tb weight, 2nd is Psat */

/* Henry's law coefficients (HLC) */
extern int nhlcdb;                              /* # number of species in the hlc database */
extern char hlcdb_chem[MXKDB][MXLFO+1];         /* formula of the species */
extern double hlcdb_dat[MXKDB][3];              /* data */
extern char hlcdb_com[MXKDB][3][MXLCOD+1];      /* comment's code for the data */

/* Hydration constants (for effective HLC) */
extern int nkhydb;                              /* # number of species in the khyd database */
extern char khydb_chem[MXKDB][MXLFO+1];         /* formula of the species */
extern double khydb_dat[MXKDB][1];              /* data */
extern char khydb_com[MXKDB][3][MXLCOD+1];      /* comment's code for the data */

/* special species list (#species) */
#define MXSPSP 3000  /* max # of special species */
extern int nspsp;                               /* # of "special" species (e.g. furane) */
extern bool lospsp[MXSPSP];                     /* logical, "true" if the species has been used. */
extern char dictsp[MXSPSP][MXLDI+1];            /* dictionary line of the special species */

/* special species chemistry (oge = Out GEnerator) */
#define MXOG 10000   /* max # of entry (reaction) in the oge database */
extern int noge;                                /* # of reactions given (oge = "out generator") */
extern int ogelab[MXOG];                        /* label for the reaction (if EXTRA or HV) */
extern char ogertve[MXOG][3][MXLFO+1];          /* reactants for reaction i */
extern char ogeprod[MXOG][MXPD][MXLFO+1];       /* products for reaction i */
extern double ogearh[MXOG][3];                  /* arrehnius coefficients (A, n, Ea) */
extern double ogestoe[MXOG][MXPD];              /* stochiometric coefficients for product j */
extern double ogeaux[MXOG][7];                  /* aux. info. for reaction i (e.g. falloff react.) */

/* fixed names data (fn database) */
#define MXFN 3000    /* max # of species with fixed names (i.e. not set by gecko) */
extern int nfn;                                 /* # of species with fixed names */
extern char namfn[MXFN][MXLCO+1];               /* names of fixed species */
extern char chemfn[MXFN][MXLFO+1];              /* standardized formula of species with fixed name */

#endif /* DATABASE_H */
