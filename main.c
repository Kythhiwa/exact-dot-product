#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <mpfr.h>


int cmp_bits(void* a, void *b, size_t size1, size_t size2) {
    char* p1 = (char*) a;
    char* p2 = (char*) b;
    int n = FLT_MANT_DIG - 1;
    int m = FLT_MANT_DIG - 1;
    if (size1 == sizeof(double)) n = DBL_MANT_DIG - 1;
    if (size1 == sizeof(long double)) n = LDBL_MANT_DIG - 1;
    if (size2 == sizeof(double)) m = DBL_MANT_DIG - 1;
    if (size2 == sizeof(long double)) m = LDBL_MANT_DIG - 1;

    int bits1[n];
    int bits2[m];
    int t = n;
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < 8; j++) {
            if (t < 0 ) break;
            bits1[--t] = ( p1[i] & (1 << j)) ? 1 : 0;
            
         }
    }
    t = m;
    for (int i = 0; i < size2; i++) {
        for (int j = 0; j < 8; j++) {
            if (t < 0) break;
            bits2[--t] = ( p2[i] & (1 << j)) ? 1 : 0;
            
         }
    }
    int f = 0;
    printf("Difference:       ");
    for (int i = 0; i < n; i++) {
        if (bits1[i] != bits2[i]) {
            f = 1;
            printf("^");
        } else {
            printf(" ");
        }
    }
    if (!f) {
        printf("\nThe mantissas are equal\n");
    } else {
        printf("\nThe mantissas are not equal\n");
    }
    return f;
}

void print_memory_binary(void *ptr, size_t size) {
    char* p = (char*) ptr;
    int n = FLT_MANT_DIG - 1;
    if (size == sizeof(double)) n = DBL_MANT_DIG - 1;
    if (size == sizeof(long double)) n = LDBL_MANT_DIG - 1;
    int bits[n];
    int t = n;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 8; j++) {
            if (t < 0 ) break;
            bits[--t] = ( p[i] & (1 << j)) ? 1 : 0;
         }
    }
    for (int i = 0; i < n; i++) printf("%d", bits[i]);
    putchar('\n');
}

int get_sign(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(x));
    return (bits >> 63) & 1;
}

int get_exp(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(x));
    return  (bits >> 52) & 0x7FF;
}

int get_gen(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(x));
    return bits & 1;
}
 
void split(double* a, double* x1, double* x2) {
    double x = *a;
    char* p = (char*) &x;
    p[0] &= 0;
    p[1] &= 0;
    p[2] &= 0;
    p[3] &= 0xE0;
    memcpy(&x, p, sizeof(double));
    *x1 = *a - x;
    *x2 = *a - *x1;
}

void insert_table(double x, double* t) {
    double c = x;
    while(t[get_exp(c) * 2 + get_gen(c)] != 0.0) {
        double e = t[get_exp(c) * 2 + get_gen(c)];
        t[get_exp(c) * 2 + get_gen(c)] = 0.0;
        c += e;

    }
    t[get_exp(c) * 2 + get_gen(c)] = c;
}

void mul(double x, double y, double* t) {
  
    double x11, x12, y1, y2;
    split(&x, &x11, &x12);
    split(&y, &y1, &y2);
    insert_table(x11 * y1, t);
    insert_table(x11 * y2, t);
    insert_table(x12 * y1, t);
    insert_table(x12 * y2, t);
}

void print_table(double* t) {
    printf("| exp|genus|  val\n");
    for (int i = 0; i <= 2048; i++) {       
        if (t[i*2] !=0) printf("|   %d|   %d|    %.52f|\n", i, 0, t[i * 2]);
        if (t[i*2+1]!=0) printf("|   %d|   %d|    %.52f|\n", i, 1, t[i * 2 + 1]);
    }
}

double two_sum_d(double* t, double a, double b) {
    double s = a + b;
    double bs = s - a;
    double as = s - bs;
    *t = (b - bs) + (a - as);
    return s;
}

double sum_table(double* t) {
    double s = 0.0, c = 0.0;
    double e;
    for (int i = 0; i <= 2048; i++) {
        for (int j = 0; j <= 1; j++) {
            if (t[i*2+j] == 0) continue;
            s = two_sum_d(&e, s, t[i * 2 + j]);
            c += e;
        }
    }
    return two_sum_d(&e, s, c) + e;
}

double random_double(double min, double max) {
    double scale = (double)rand() / RAND_MAX;
    return min + scale * (max - min);
}

void generate(double* x, double* y, int n, double l, double r) {
    for (int i = 0; i < n; i++) {
        x[i] = random_double(l, r);
        y[i] = random_double(l, r);
    }
}

double mpfr(double* x, double* y, int n) {
    mpfr_t exact_result;
    mpfr_init2(exact_result, 256);

    mpfr_set_d(exact_result, 0.0, MPFR_RNDN);
    for (int i = 0; i < n; i++) {
        mpfr_t mpfr_x, term;
        
        mpfr_init2(mpfr_x, 256);
        mpfr_set_d(mpfr_x, x[i], MPFR_RNDN); 
        
        mpfr_init2(term, 256);
        mpfr_mul_d(term, mpfr_x, y[i], MPFR_RNDN);
        
        
        mpfr_add(exact_result, exact_result, term, MPFR_RNDN);
        
        
        mpfr_clears(mpfr_x, term, NULL);    
    }

    double exact_double = mpfr_get_d(exact_result, MPFR_RNDN);
    return exact_double;
}

void shuffle(double* x, double* y, int n) {
    for (int i = n - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        double tx = x[i], ty = y[i];
        x[i] = x[j]; y[i] = y[j];
        x[j] = tx; y[j] = ty;
    }
}


int test(double* x, double* y, int n) {
    double* t = (double*) calloc(2049 * 2, sizeof(double));
    double* to = (double*) calloc(2049 * 2, sizeof(double));
    double* tO = (double*) calloc(2049 * 2, sizeof(double));
    long double exact = 0L;

    int k = 1024;
    for (int i = 0; i < n; i++) {
        int exp_x = get_exp(x[i]) - 1023;
        int exp_y = get_exp(y[i]) - 1023;
        int max_exp = (exp_x > exp_y) ? exp_x : exp_y;
        int min_exp = (exp_x < exp_y) ? exp_y : exp_x;
        if (max_exp > 485) {
            if (exp_x > exp_y) {
                mul(ldexp(x[i], -k), y[i], to);
            } else {
                mul(x[i], ldexp(y[i], -k), to);
            }
        } else if (min_exp < -485) {
            if (exp_x < exp_y) {
                mul(ldexp(x[i], k), y[i], tO);
            } else {
                mul(x[i], ldexp(y[i], k), tO);
            }
        } else {
            mul(x[i], y[i], t);
        }
        long double w =  (long double) x[i] *  (long double) y[i];
        double x1, x2;
        exact += w;

    }
    double s_exact = (double) exact;
    double s = sum_table(t);
    double so = ldexp(sum_table(to), k);
    double sO = ldexp(sum_table(tO), -k);
    double c;
    double p = 0.0;
    s += two_sum_d(&c, s, so);
    p += c;
    s += two_sum_d(&c, s, sO);
    p += c;
    s += p;
    double* t1 = (double*) calloc(2049 * 2, sizeof(double));
    double* to1 = (double*) calloc(2049 * 2, sizeof(double));
    double* tO1 = (double*) calloc(2049 * 2, sizeof(double));

    shuffle(x, y, n);
    for (int i = 0; i < n; i++) {
        int exp_x = get_exp(x[i]) - 1023;
        int exp_y = get_exp(y[i]) - 1023;
        int max_exp = (exp_x > exp_y) ? exp_x : exp_y;
        int min_exp = (exp_x < exp_y) ? exp_y : exp_x;
        if (max_exp > 485) {
            if (exp_x > exp_y) {
                mul(ldexp(x[i], -k), y[i], to1);
            } else {
                mul(x[i], ldexp(y[i], -k), to1);
            }
        } else if (min_exp < -485) {
            if (exp_x < exp_y) {
                mul(ldexp(x[i], k), y[i], tO1);
            } else {
                mul(x[i], ldexp(y[i], k), tO1);
            }
        } else {
            mul(x[i], y[i], t1);
        }
    }
    double s1 = sum_table(t1);
    double so1 = ldexp(sum_table(to1), k);
    double sO1 = ldexp(sum_table(tO1), -k);
    s1 += two_sum_d(&c, so1, sO1);
    s1 += c;

    double ss = mpfr(x, y, n);
    printf("Long double:      ");
    print_memory_binary(&s_exact, sizeof(s_exact));
    printf("Computed:         ");
    print_memory_binary(&s, sizeof(s));    
    printf("Exact mpfr:       ");
    print_memory_binary(&ss, sizeof(ss));
    int f1 = cmp_bits(&s, &ss, sizeof(s), sizeof(ss));

    printf("Computed:         ");
    print_memory_binary(&s, sizeof(s));
    printf("Shuffled:         ");
    print_memory_binary(&s1, sizeof(s1));
    f1 += cmp_bits(&s, &s1, sizeof(s), sizeof(s1));
    free(t);
    free(to);
    free(tO);
    free(t1);
    free(to1);
    free(tO1);
    return f1 != 0;
}

int main() {
    srand(7);
    int n = 100;
    double x[n], y[n];
    int cnt = 0;
    int t1 = 1000;
    for (int  i = 0; i < t1;i++) {
        printf("\nTest %d\n", i);
        generate(x, y, n, -1000, 1000);
        
        printf("\n");
        int f2 = test(x, y, n);
        if (f2) cnt++;
    }
    int t2 = 1000;
    int cnt1 = 0;
    for (int  i = 0; i < t2; i++) {
        printf("\nTest %d\n", i);
        generate(x, y, n, -1e-152, 1e-152);
        int f = test(x, y, n);
        if (f) cnt1++;
    }  
    int t3 = 1000;
    int cnt2 = 0;
    for (int  i = 0; i < t3; i++) {
        printf("\nTest %d\n", i);
        generate(x, y, n, -1e152, 1e152);
        int f = test(x, y, n);
        if (f) cnt2++;
    }  
    printf("%d/%d\n", t1 - cnt, t1);
    printf("%d/%d\n", t2 - cnt1, t2);
    printf("%d/%d\n", t3 - cnt2, t2);
}
