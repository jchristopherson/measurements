#ifndef C_ARRAY_H_
#define C_ARRAY_H_

#ifdef __cplusplus
extern "C" {
#endif

void create_random_array(int n, double *x);
double array_sum(int n, const double *x);
void array_sort(int n, double *x);


#ifdef __cplusplus
}
#endif
#endif
