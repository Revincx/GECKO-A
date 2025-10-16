/* Key flags and parameters for GECKO-A */
#include "keyflag.h"
#include "keyparameter.h"
#include <string.h>

/* Global variable definitions */
double critvp;
double brcut;
double yldcut;
double rxloss;
int maxgen;
double TK;

/* SAR selectors */
int kisom_sar;
int kdiss_sar;
int kohadd_sar;
int pvap_sar;

/* FLAGS */
bool g2pfg;
bool g2wfg;
bool isomerfg;
bool enolflg;
bool pamfg;
bool highnoxfg;
bool rx_ro_no2;
double C_NO2;
bool rx_ro2_no2;
bool dhffg;
bool chafg;
bool bimolecrx4criegee;
bool rx_ro2_multiclass;
bool rx_ro2_oh;

/* INFOS */
bool sar_only_fg;
bool wrtsarinfo;

/* OUTPUT FLAGS */
bool screenfg;
bool wrtdhf;
bool wrtpvap;
bool wrthenry;
bool wrtref;
bool wrttg;
bool wrtdepo;
bool wrtkivoci;
bool wrtmaxyield;

/* Set default values for all parameters */
void define_defaults(void) {
    /* various cut-offs */
    maxgen = 20;              /* maximum # of generations allowed */
    TK     = 298.;            /* temperature for which to generate the mechanism */
    critvp = -13;             /* log(Pvap) below which chemistry is ignored */
    brcut  = 0.05;            /* rx cut off branching ratio below which a reaction pathway is ignored */
    yldcut = 1E-3;            /* rx cut off if: (rx branching ratio)*(yield) is below threshold */
    rxloss = 1E-10;           /* send to lost carbon if species yield below threshold */
    C_NO2  = 2.5E10;          /* NO2 concentration (molec.cm-3) used in roselector (rochem) for RO+NO2 reaction */

    /* reaction switches */
    rx_ro_no2 = false;        /* add the RO+NO2 reaction to the RO chemistry */
    rx_ro2_no2 = false;       /* add the RO2+NO2 reversible reaction to the RO2 chemistry */
    rx_ro2_multiclass = true; /* treat RO2+RO2 with the 9 classes (true) of ro2 or only CH3O2 (false) */

    /* SAR */
    pvap_sar   = 2;           /* SAR for vapor pressure: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1 */
    kisom_sar  = 2;           /* R(O.) selector for isomerisation reaction rate (1=Atkinson 2007; 2=Vereecken 2009) */
    kdiss_sar  = 2;           /* R(O.) selector for decomposition reaction rate (1=Atkinson 2007; 2=Vereecken 2009) */
    kohadd_sar = 3;           /* OH addition to alkene (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016) */
    
    /* flags */
    g2pfg      = true;        /* flag to write for gas/particle partitioning "reactions" */
    g2wfg      = false;       /* flag to write for gas/wall partitioning "reactions" (chamber simulations) */
    isomerfg   = true;        /* allow isomer substitution */
    highnoxfg  = false;       /* only NO chemistry for RO2 considered */
    dhffg      = false;       /* flag to activate DHF formation */
    chafg      = false;       /* consider cyclic hemi acetal (CHA) species as non volatile */
    pamfg      = false;       /* "PAM" simulations: H chemistry could be considered */
    bimolecrx4criegee = false; /* allow bimolecular reaction for stabilized criegee (keep H2O) */
    enolflg    = true;        /* allow to switch enols to ketones */
    rx_ro2_oh  = false;       /* add the RO2+OH reaction in RO2 chemistry */

    /* Debug/infos */
    wrtsarinfo = false;       /* write info about SAR (groups, ...) in dedicated files */
    sar_only_fg = false;      /* write kOH, kNO3, kO3 for SAR assessment - if true, run the reactions for the input species only */

    /* default directories */
    strcpy(dirgecko, "../");
    strcpy(dirout, "OUT/");

    /* write output files */
    screenfg   = false;       /* write additional info on screen during generation */
    wrtdhf     = false;       /* write output file with formation enthalpies information */
    wrtpvap    = true;        /* write output file with the vapor pressures */
    wrthenry   = true;        /* write output file with the Henry's law coefs */
    wrtref     = true;        /* write the references for the reactions in the mechanism */
    wrttg      = true;        /* write output file with Tg data */
    wrtdepo    = false;       /* write output file with data for deposition */
    wrtkivoci  = false;       /* write output file with kivoci and kjvocj */
    wrtmaxyield = false;      /* write output file with maximum yields */
}
