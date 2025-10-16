/* Mini dictionary module for GECKO-A */
#ifndef MINIDICT_H
#define MINIDICT_H

#include "keyparameter.h"

#define MXMINID 3000

/* Global data */
extern int nmini;
extern char nam_minid[MXMINID][MXLCO+1];
extern char fo_minid[MXMINID][MXLFO+1];

/* Function declarations */
void get_fo(const char *nam, char *formula);
void add_fo(const char *nam, const char *formula);
void clean_minid(void);

#endif /* MINIDICT_H */
