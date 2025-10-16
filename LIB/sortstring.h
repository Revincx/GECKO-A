/* String sorting module using quicksort + insertion sort */
#ifndef SORTSTRING_H
#define SORTSTRING_H

#define QSORT_THRESHOLD 16

/* Sort an array of strings in alphabetical order */
void sort_string(char **strings, int count, int max_len);

/* Sort a 2D array of strings */
void sort_string_array(char strings[][256], int count);

#endif /* SORTSTRING_H */
