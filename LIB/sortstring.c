/* String sorting module using quicksort + insertion sort 
 * Based on optimized combination of QSORT and INSERTION sorting
 */
#include "sortstring.h"
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

#define QSORT_STACK_SIZE 32

typedef struct {
    int low;
    int high;
} qsort_stack_t;

/* Helper: compare two strings */
static inline bool less_than(const char *a, const char *b) {
    return strcmp(a, b) < 0;
}

/* Helper: swap two strings */
static void swap_strings(char *a, char *b, int max_len) {
    char *temp = (char *)malloc(max_len + 1);
    strcpy(temp, a);
    strcpy(a, b);
    strcpy(b, temp);
    free(temp);
}

/* Helper: circular shift from left to right */
static void r_shift(char **strings, int left, int right, int max_len) {
    char *hold = (char *)malloc(max_len + 1);
    strcpy(hold, strings[right]);
    
    for (int i = right; i > left; i--) {
        strcpy(strings[i], strings[i-1]);
    }
    strcpy(strings[left], hold);
    free(hold);
}

/* Sort an array of strings */
void sort_string(char **strings, int count, int max_len) {
    if (count <= 1) return;
    
    qsort_stack_t stack[QSORT_STACK_SIZE];
    int stack_top = 0;
    int low, high, left, right, mid;
    int right_size, left_size;
    
    /* Quick sort for large partitions */
    if (count > QSORT_THRESHOLD) {
        low = 0;
        high = count - 1;
        
        while (true) {
            mid = (low + high) / 2;
            
            if (less_than(strings[mid], strings[low])) {
                swap_strings(strings[mid], strings[low], max_len);
            }
            if (less_than(strings[high], strings[mid])) {
                swap_strings(strings[high], strings[mid], max_len);
                if (less_than(strings[mid], strings[low])) {
                    swap_strings(strings[mid], strings[low], max_len);
                }
            }
            
            left = low + 1;
            right = high - 1;
            
            /* Collapse walls */
            while (true) {
                while (less_than(strings[left], strings[mid])) {
                    left++;
                }
                while (less_than(strings[mid], strings[right])) {
                    right--;
                }
                
                if (left < right) {
                    swap_strings(strings[left], strings[right], max_len);
                    if (mid == left) {
                        mid = right;
                    } else if (mid == right) {
                        mid = left;
                    }
                    left++;
                    right--;
                } else {
                    if (left == right) {
                        left++;
                        right--;
                    }
                    break;
                }
            }
            
            /* Set up indices for next iteration */
            right_size = right - low;
            left_size = high - left;
            
            if (right_size <= QSORT_THRESHOLD) {
                if (left_size <= QSORT_THRESHOLD) {
                    if (stack_top < 1) break;
                    stack_top--;
                    low = stack[stack_top].low;
                    high = stack[stack_top].high;
                } else {
                    low = left;
                }
            } else if (left_size <= QSORT_THRESHOLD) {
                high = right;
            } else if (right_size > left_size) {
                stack[stack_top].low = low;
                stack[stack_top].high = right;
                stack_top++;
                low = left;
            } else {
                stack[stack_top].low = left;
                stack[stack_top].high = high;
                stack_top++;
                high = right;
            }
        }
    }
    
    /* Insertion sort for remaining small partitions */
    low = 0;
    high = count - 1;
    
    /* Find smallest in first QSORT_THRESHOLD elements */
    left = low;
    int limit = (low + QSORT_THRESHOLD < high) ? low + QSORT_THRESHOLD : high;
    for (right = low + 1; right <= limit; right++) {
        if (less_than(strings[right], strings[left])) {
            left = right;
        }
    }
    if (left != low) {
        swap_strings(strings[left], strings[low], max_len);
    }
    
    /* Insertion sort from left to right */
    for (right = low + 2; right <= high; right++) {
        left = right - 1;
        if (less_than(strings[right], strings[left])) {
            while (less_than(strings[right], strings[left-1])) {
                left--;
            }
            r_shift(strings, left, right, max_len);
        }
    }
}

/* Sort a 2D array of strings (convenience wrapper) */
void sort_string_array(char strings[][256], int count) {
    if (count <= 1) return;
    
    /* Create array of pointers for easier manipulation */
    char **ptrs = (char **)malloc(count * sizeof(char *));
    for (int i = 0; i < count; i++) {
        ptrs[i] = strings[i];
    }
    
    sort_string(ptrs, count, 256);
    free(ptrs);
}
