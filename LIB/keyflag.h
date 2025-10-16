/* Key flags and parameters for GECKO-A */
#ifndef KEYFLAG_H
#define KEYFLAG_H

#include <stdbool.h>

/* various cut-offs */
extern double critvp;     /* log(Pvap) below which chemistry is ignored */
extern double brcut;      /* rx cut off branching ratio below which a reaction pathway is ignored */
extern double yldcut;     /* rx cut off if: (rx branching ratio)*(yield) is below threshold */
extern double rxloss;     /* send to lost carbon if species yield below threshold */

extern int maxgen;        /* maximum # of generations allowed */
extern double TK;         /* user-set temperature */

/* reference values */
#define TREF 298.0        /* reference temperature */
#define MREF 2.5E19       /* default 3rd body M (molec cm-3) */

/* CHEMISTRY FLAG - SAR SELECTOR */
extern int kisom_sar;     /* R(O.) selector for isomerisation reaction rate (1=Atkinson 2007; 2=Vereecken 2009) */
extern int kdiss_sar;     /* R(O.) selector for decomposition reaction rate (1=Atkinson 2007; 2=Vereecken 2009) */
extern int kohadd_sar;    /* OH addition to alkene (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016) */
extern int pvap_sar;      /* select the SAR for vapor pressure: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1 */

/* FLAGS */
extern bool g2pfg;                    /* flag to write for gas/particle partitioning "reactions" */
extern bool g2wfg;                    /* flag to write for gas/wall partitioning "reactions" (chamber simulations) */
extern bool isomerfg;                 /* allow isomer substitution (O=no substitution, 1=allow) */
extern bool enolflg;                  /* allow to switch enols to ketones */
extern bool pamfg;                    /* "PAM" simulations: H chemistry could be considered */
extern bool highnoxfg;                /* only NO chemistry for RO2 considered */
extern bool rx_ro_no2;                /* add the RO+NO2 reaction to the RO chemistry */
extern double C_NO2;                  /* NO2 concentration (molec.cm-3) used in roselector (rochem) for RO+NO2 reaction */
extern bool rx_ro2_no2;               /* add the RO2+NO2 reversible reaction to the RO2 chemistry */
extern bool dhffg;                    /* flag to activate DHF formation */
extern bool chafg;                    /* consider cyclic hemi acetal (CHA) species as non volatile */
extern bool bimolecrx4criegee;        /* allow bimolecular reaction for stabilized criegee (keep H2O) */
extern bool rx_ro2_multiclass;        /* treat RO2+RO2 with the 9 classes (true) of ro2 or only CH3O2 (false) */
extern bool rx_ro2_oh;                /* add the RO2+OH reaction in RO2 chemistry */

/* INFOS */
extern bool sar_only_fg;              /* write kOH, kNO3, kO3 for SAR assessment - if true, run the reactions for the input species only */
extern bool wrtsarinfo;               /* write info about SAR (groups, ...) in dedicated files */
#define WRTOPEINFO false              /* write "operator info" in the output files - OLD STUFF */

/* OUTPUT FLAGS */
extern bool screenfg;                 /* write additional info on screen during generation */
extern bool wrtdhf;                   /* write output file with formation enthalpies information */
extern bool wrtpvap;                  /* write output file with the vapor pressures */
extern bool wrthenry;                 /* write output file with the Henry's law coefs */
extern bool wrtref;                   /* write the references for the reactions in the mechanism */
extern bool wrttg;                    /* write output file with Tg data */
extern bool wrtdepo;                  /* write output file with data for deposition */
extern bool wrtkivoci;                /* write output file with kivoci and kjvocj */
extern bool wrtmaxyield;              /* write output file with maximum yields */

/* Function to set default values */
void define_defaults(void);

#endif /* KEYFLAG_H */
