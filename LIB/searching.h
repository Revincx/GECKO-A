/* Searching utilities for GECKO-A */
#ifndef SEARCHING_H
#define SEARCHING_H

/* Binary tree search for string
 * Returns: <0 if not found (|result| == pointer where aseek should be inserted)
 *          >0 if found (|result| == pointer where aseek is inserted in alist)
 */
int srh5(const char *aseek, char **alist, int nlist);

/* Binary tree search for formula in dictionary
 * Returns: <0 if not found (|result| == pointer where chem should be inserted)
 *          >0 if found (|result| == pointer where chem is inserted in dict)
 */
int srch(int nrec, const char *chem, char **dict);

/* Return the rank (position) of an integer in a sorted list of integers
 * (from large to small number)
 */
int search_ipos(int iseek, const int *ilist, int list_size);

#endif /* SEARCHING_H */
