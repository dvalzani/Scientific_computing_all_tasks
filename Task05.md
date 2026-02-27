### Task05

## Defined the following vector: vec = [1.0, 1.0e16, -1.0e16, -0.5]

Code for calculating the sum of this vector in 3 different ways:

- loop

- gsl_vector_sum

- kahan algorithm

  ```c
  #include <stdio.h>
  #include <math.h>
  #include <gsl/gsl_vector.h>

  // Kahan summation
  double kahan_sum(const double *arr, size_t n) {
    double sum = 0.0;
    double c = 0.0;  // compensazione
    for (size_t i = 0; i < n; ++i) {
        double y = arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
  }

  int main(void) {
    // Vettore dell'esercizio
    double vec[4] = {1.0, 1.0e16, -1.0e16, -0.5};
    const size_t N = 4;
    const double true_sum = 0.5;

    // 1) Somma con for loop semplice
    double sum_loop = 0.0;
    for (size_t i = 0; i < N; ++i) {
        sum_loop += vec[i];
    }

  // 2) Somma con funzione di libreria GSL (modo corretto):
  // GSL NON ha "gsl_vector_sum", quindi facciamo una somma manuale
  gsl_vector_const_view vview = gsl_vector_const_view_array(vec, N);
  const gsl_vector *v = &vview.vector;

  double sum_gsl = 0.0;
  for (size_t i = 0; i < N; ++i) {
    sum_gsl += gsl_vector_get(v, i);
  }

    // 3) Somma con Kahan
    double sum_kahan = kahan_sum(vec, N);

    // Stampa risultati
    printf("For loop sum      = %.16f\n", sum_loop);
    printf("GSL gsl_vector_sum= %.16f\n", sum_gsl);
    printf("Kahan sum         = %.16f\n", sum_kahan);
    printf("True (analytical) = %.16f\n", true_sum);

    // Errori relativi
    double rel_err_loop  = fabs(sum_loop  - true_sum) / fabs(true_sum);
    double rel_err_gsl   = fabs(sum_gsl   - true_sum) / fabs(true_sum);
    double rel_err_kahan = fabs(sum_kahan - true_sum) / fabs(true_sum);

    printf("\nRelative error (Loop) : %.2e\n", rel_err_loop);
    printf("Relative error (GSL)  : %.2e\n", rel_err_gsl);
    printf("Relative error (Kahan): %.2e\n", rel_err_kahan);

    return 0;
  }

  ```

 ## The output is:
```
For loop sum      = -0.5000000000000000
GSL gsl_vector_sum= -0.5000000000000000
Kahan sum         = -0.5000000000000000
True (analytical) = 0.5000000000000000
  
Relative error (Loop) : 2.00e+00
Relative error (GSL)  : 2.00e+00
Relative error (Kahan): 2.00e+00
```

### Task 05 (a) – Are the three results the same? Why / why not?

The three results calculated with the 3 different methods are the same, but they are different from the real value.

This happens because the "1 value" added to 10^16, is lost inside the precision of the double number

# One possible solution is:

- Use longdouble instead double:

```
  Loop (long double)  = 0.50000000000000000000
  Kahan (long double) = 0.50000000000000000000
  True                = 0.50000000000000000000
```

## B) code for daxpy alghorithm with x and y vectors filled with a gaussian random variables with mean 0 and std 1

```c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

int main(void) {
    const size_t N = 1000000;    // dimensione vettori
    const double a = 3.0;        // stesso 'a' delle altre task
    const unsigned long seed = 123; // seed per riproducibilità

    // 1) Alloco i vettori
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));

    if (!x || !y || !d) {
        fprintf(stderr, "Errore: malloc fallita.\n");
        free(x); free(y); free(d);
        return 1;
    }

    // 2) Inizializzo il generatore di numeri casuali GSL
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    if (!r) {
        fprintf(stderr, "Errore: gsl_rng_alloc fallita.\n");
        free(x); free(y); free(d);
        return 1;
    }
    gsl_rng_set(r, seed);

    // 3) Genero x, y ~ N(0,1) e calcolo d = a*x + y
    for (size_t i = 0; i < N; ++i) {
        x[i] = gsl_ran_gaussian(r, 1.0);  // mean=0, std=1
        y[i] = gsl_ran_gaussian(r, 1.0);
        d[i] = a * x[i] + y[i];
    }

    // -------------------------------
    // TEST 1: Verifica numerica element-wise
    // -------------------------------
    double max_abs_err = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double ref = a * x[i] + y[i];        // "teorico"
        double err = fabs(d[i] - ref);
        if (err > max_abs_err) {
            max_abs_err = err;
        }
    }
    printf("Max absolute difference between d and a*x + y: %.3e\n", max_abs_err);

    // -------------------------------
    // TEST 2: Verifica statistica
    // -------------------------------
    // Se x,y ~ N(0,1) indipendenti:
    // E[x] = E[y] = 0; Var(x)=Var(y)=1
    // d = a*x + y -> E[d]=0, Var(d) = a^2 + 1

    double mean_x = gsl_stats_mean(x, 1, N);
    double std_x  = gsl_stats_sd(x, 1, N);

    double mean_y = gsl_stats_mean(y, 1, N);
    double std_y  = gsl_stats_sd(y, 1, N);

    double mean_d = gsl_stats_mean(d, 1, N);
    double std_d  = gsl_stats_sd(d, 1, N);

    double std_d_theoretical = sqrt(a*a + 1.0);

    printf("\n--- Statistics ---\n");
    printf("mean(x) = %.3e, std(x) = %.3e  (expected: 0, 1)\n",
           mean_x, std_x);
    printf("mean(y) = %.3e, std(y) = %.3e  (expected: 0, 1)\n",
           mean_y, std_y);
    printf("mean(d) = %.3e\n", mean_d);
    printf("std(d)  = %.3e, expected std(d) = %.3e\n",
           std_d, std_d_theoretical);

    // Piccolo check "booleano"
    double tol_mean = 5e-3;
    double tol_std  = 5e-3;

    int ok_mean = fabs(mean_d - 0.0) < tol_mean;
    int ok_std  = fabs(std_d - std_d_theoretical) < tol_std;

    printf("\n--- Checks ---\n");
    printf("|mean(d) - 0| < %.1e ? -> %s\n",
           tol_mean, ok_mean ? "TRUE" : "FALSE");
    printf("|std(d) - sqrt(a^2+1)| < %.1e ? -> %s\n",
           tol_std, ok_std ? "TRUE" : "FALSE");

    // Pulizia
    gsl_rng_free(r);
    free(x);
    free(y);
    free(d);

    return 0;
}

```

# The Output
```
Max absolute difference between d and a*x + y: 0.000e+00

--- Statistics ---
mean(x) = 4.537e-04, std(x) = 1.000e+00  (expected: 0, 1)
mean(y) = -1.504e-03, std(y) = 9.993e-01  (expected: 0, 1)
mean(d) = -1.426e-04
std(d)  = 3.162e+00, expected std(d) = 3.162e+00

--- Checks ---
|mean(d) - 0| < 5.0e-03 ? -> TRUE
|std(d) - sqrt(a^2+1)| < 5.0e-03 ? -> TRUE

```

**Makefile for both the codes**

```c
CC      := gcc
CFLAGS  := -O2 -Wall -Wextra -std=c11
LDLIBS  := -lgsl -lgslcblas -lm

TARGETS := a_floating_point b_daxpy_gaussian

all: $(TARGETS)

a_floating_point: a_floating_point.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

b_daxpy_gaussian: b_daxpy_gaussian.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: all clean


```



  
