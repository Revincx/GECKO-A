/* Log writing tools for GECKO-A */
#include "logtool.h"
#include "keyparameter.h"
#include "keyflag.h"
#include <stdio.h>

/* Write operating conditions and flags to log file */
void wrtlog(void) {
    FILE *logfile;
    char filename[256];
    
    /* Open log file - using LOGU file unit number would map to a specific file */
    sprintf(filename, "%slog.txt", dirout);
    logfile = fopen(filename, "a");
    if (!logfile) {
        fprintf(stderr, "Warning: Could not open log file %s\n", filename);
        logfile = stdout;  /* Fall back to stdout */
    }
    
    /* Write mechanism parameters to log */
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, " GECKO VERSION=v0.0.0 \n");
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, " maxgen = %d\n", maxgen);
    fprintf(logfile, "  mechanism temperature (K) used for branching ratio: %6.1f\n", TK);
    fprintf(logfile, "  reference temperature (K): %6.1f\n", TREF);
    fprintf(logfile, "  default 3rd body M (molec cm-3): %10.2e\n", MREF);
    fprintf(logfile, "          \n");
    
    /* Mechanism options */
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, "    ----  MECHANISM  ---- \n");
    fprintf(logfile, " ---------------------------- \n");
    
    if (g2pfg) {
        fprintf(logfile, " gas/particle mass transfer: included\n");
    } else {
        fprintf(logfile, " gas/particle mass transfer: not included\n");
    }
    
    if (g2wfg) {
        fprintf(logfile, " gas/wall mass transfer: included\n");
    } else {
        fprintf(logfile, " gas/wall mass transfer: not included\n");
    }
    
    if (dhffg) {
        fprintf(logfile, " DHF isomerisation: yes\n");
    } else {
        fprintf(logfile, " DHF isomerisation: no\n");
    }
    fprintf(logfile, "          \n");
    
    /* SAR options */
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, "    --------- SAR -------- \n");
    fprintf(logfile, " ---------------------------- \n");
    
    if (pvap_sar == 1) {
        fprintf(logfile, " vapor pressure scheme=M&Y \n");
    } else if (pvap_sar == 2) {
        fprintf(logfile, " vapor pressure scheme=Nannoolal 2008\n");
    } else if (pvap_sar == 3) {
        fprintf(logfile, " vapor pressure scheme=SIMPOL-1 \n");
    }
    
    if (kohadd_sar == 1) {
        fprintf(logfile, " OH addition SAR =Peeters 1997 \n");
    } else if (kohadd_sar == 2) {
        fprintf(logfile, " OH addition SAR =Ziemann 2009 \n");
    } else if (kohadd_sar == 3) {
        fprintf(logfile, " OH addition SAR =Jenkin 2018  \n");
    }
    
    if (kisom_sar == 1) {
        fprintf(logfile, " Alkoxy isomerisation: Atkinson 2007 \n");
    } else if (kisom_sar == 2) {
        fprintf(logfile, " Alkoxy isomerisation: Vereecken 2009 \n");
    }
    
    if (kdiss_sar == 1) {
        fprintf(logfile, " Alkoxy dissociation: Atkinson 2007 \n");
    } else if (kdiss_sar == 2) {
        fprintf(logfile, " Alkoxy dissociation: Vereecken 2009 \n");
    }
    
    fprintf(logfile, "          \n");
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, "   -- OPTIONAL CHEMISTRY --   \n");
    fprintf(logfile, " ---------------------------- \n");
    
    if (rx_ro2_oh) {
        fprintf(logfile, " RO2 + OH reactions activated\n");
        fprintf(logfile, "    >> Warning : chemistry for (OOOH) is simple here and needs to be improved\n");
    }
    
    if (rx_ro2_no2) {
        fprintf(logfile, " RO2 + NO2 reactions activated\n");
        fprintf(logfile, "    >> Warning : species bearing (OONO2) only react by ROONO2 -> RO2 + NO2 decomposition\n");
    }
    
    if (rx_ro_no2) {
        fprintf(logfile, " RO + NO2 -> R(ONO2) reactions activated\n");
        fprintf(logfile, "    >> Warning : rate calculated with [NO2]=%10.2e molec.cm-3 to eliminate minor pathways of RO reactions\n", C_NO2);
    }
    
    if (bimolecrx4criegee) {
        fprintf(logfile, " Bimolecular reaction for stabilized criegee allowed (keep H2O)\n");
    }
    fprintf(logfile, "          \n");
    
    /* Reductions */
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, "    ----- Reductions ----- \n");
    fprintf(logfile, " ---------------------------- \n");
    fprintf(logfile, " critical vapor pressure (atm) = %f\n", critvp);
    
    if (isomerfg) {
        fprintf(logfile, " isomerisation allowed \n");
    } else {
        fprintf(logfile, " no isomerisation \n");
    }
    
    fprintf(logfile, " Cut off branching ratio below which a reaction pathway is ignored: %f\n", brcut);
    fprintf(logfile, " high-NOx flag: %d\n", highnoxfg ? 1 : 0);
    fprintf(logfile, " All classes of RO2 in RO2+RO2: %d\n", rx_ro2_multiclass ? 1 : 0);
    fprintf(logfile, "          \n");
    fprintf(logfile, " ---------------------------- \n");
    
    if (logfile != stdout) {
        fclose(logfile);
    }
}
