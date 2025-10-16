/* Database module for GECKO-A */
#include "database.h"

/* VOC+OH database */
int nkohdb = 0;
char kohdb_chem[MXKDB][MXLFO+1];
double kohdb_298[MXKDB];
double kohdb_arr[MXKDB][3];
char kohdb_com[MXKDB][3][MXLCOD+1];

/* VOC+O3 database */
int nko3db = 0;
char ko3db_chem[MXKDB][MXLFO+1];
double ko3db_298[MXKDB];
double ko3db_arr[MXKDB][3];
char ko3db_com[MXKDB][3][MXLCOD+1];

/* VOC+NO3 database */
int nkno3db = 0;
char kno3db_chem[MXKDB][MXLFO+1];
double kno3db_298[MXKDB];
double kno3db_arr[MXKDB][3];
char kno3db_com[MXKDB][3][MXLCOD+1];

/* VOC+OH mechanism */
int nkwoh = 0;
char kwoh_rct[MXKWR][MXLFO+1];
int nkwoh_pd[MXKWR];
double kwoh_yld[MXKWR][MXNR];
char kwoh_pd[MXKWR][MXNR][MXLFO+1];
char kwoh_copd[MXKWR][MXNR][MXCOPD][MXLCO+1];
char kwoh_com[MXKWR][3][MXLCOD+1];

/* VOC+O3 mechanism */
int nkwo3 = 0;
char kwo3_rct[MXKWR][MXLFO+1];
int nkwo3_pd[MXKWR];
double kwo3_yld[MXKWR][MXNR];
char kwo3_pd[MXKWR][MXNR][MXLFO+1];
char kwo3_copd[MXKWR][MXNR][MXCOPD][MXLCO+1];
char kwo3_com[MXKWR][3][MXLCOD+1];

/* VOC+NO3 mechanism */
int nkwno3 = 0;
char kwno3_rct[MXKWR][MXLFO+1];
int nkwno3_pd[MXKWR];
double kwno3_yld[MXKWR][MXNR];
char kwno3_pd[MXKWR][MXNR][MXLFO+1];
char kwno3_copd[MXKWR][MXNR][MXCOPD][MXLCO+1];
char kwno3_com[MXKWR][3][MXLCOD+1];

/* RO2 chemistry */
int nkwro2 = 0;
double kwro2_arrh[MXKWR][3];
double kwro2_stoi[MXKWR][4];
char kwro2_rct[MXKWR][2][MXLFO+1];
char kwro2_prd[MXKWR][4][MXLFO+1];
char kwro2_com[MXKWR][4][MXLCOD+1];

/* RCOO2 chemistry */
int nkwrco3 = 0;
double kwrco3_arrh[MXKWR][3];
double kwrco3_stoi[MXKWR][4];
char kwrco3_rct[MXKWR][2][MXLFO+1];
char kwrco3_prd[MXKWR][4][MXLFO+1];
char kwrco3_com[MXKWR][4][MXLCOD+1];

/* RO chemistry */
int nkwro = 0;
double kwro_arrh[MXKWR][3];
double kwro_stoi[MXKWR][4];
char kwro_rct[MXKWR][2][MXLFO+1];
char kwro_prd[MXKWR][4][MXLFO+1];
char kwro_com[MXKWR][4][MXLCOD+1];

/* Criegee chemistry */
int nkwcri = 0;
double kwcri_arrh[MXKWR][3];
double kwcri_stoi[MXKWR][4];
char kwcri_rct[MXKWR][2][MXLFO+1];
char kwcri_prd[MXKWR][4][MXLFO+1];
char kwcri_com[MXKWR][3][MXLCOD+1];

/* photolysis labels & known photolysis reactions */
int njdat = 0;
char jchem[MXKWR][MXLFO+1];
char jprod[MXKWR][2][MXLFO+1];
int jlabel[MXKWR];
char coprodj[MXKWR][MXLCO+1];
char j_com[MXKWR][3][MXLCOD+1];
int nj40 = 0;
int jlab40[MXKWR];
double j40[MXKWR];

/* Benson groups database */
int nbson = 0;
double bsonval[MXBG];
char bsongrp[MXBG][LBG+1];

/* Simpol SAR for vapor pressure */
double p0simpgrp298 = 0.0;
double psimpgrp298[30];
double h0simpgrp298 = 0.0;
double hsimpgrp298[30];

/* Nannoolal SAR for vapor pressure and Tb */
double nanweight[219][2];

/* Henry's law coefficients (HLC) */
int nhlcdb = 0;
char hlcdb_chem[MXKDB][MXLFO+1];
double hlcdb_dat[MXKDB][3];
char hlcdb_com[MXKDB][3][MXLCOD+1];

/* Hydration constants (for effective HLC) */
int nkhydb = 0;
char khydb_chem[MXKDB][MXLFO+1];
double khydb_dat[MXKDB][1];
char khydb_com[MXKDB][3][MXLCOD+1];

/* special species list (#species) */
int nspsp = 0;
bool lospsp[MXSPSP];
char dictsp[MXSPSP][MXLDI+1];

/* special species chemistry (oge = Out GEnerator) */
int noge = 0;
int ogelab[MXOG];
char ogertve[MXOG][3][MXLFO+1];
char ogeprod[MXOG][MXPD][MXLFO+1];
double ogearh[MXOG][3];
double ogestoe[MXOG][MXPD];
double ogeaux[MXOG][7];

/* fixed names data (fn database) */
int nfn = 0;
char namfn[MXFN][MXLCO+1];
char chemfn[MXFN][MXLFO+1];
