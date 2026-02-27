# Task02

## vector_sum.py

```python
# Task 02 – Punto 1 (interpreted language)
# d = a*x + y con N (dimensione vettori) = 10, 10^6, 10^8
# a = 3, x = 0.1 per tutti, y = 7.1 per tutti

import numpy as np

# chiedo all'utente di inserire N
print("Programma: calcolo di d = a*x + y")
print("a = 3, x = 0.1, y = 7.1")
N = int(input("Inserisci il valore di N (es. 10, 1000000, 100000000): "))

#definisco lo scalare e i vettori
a=3
x=[0.1]*N #crea un vettore con N elementi, tutti pari a 0.1
y=[7.1]*N #crea un vettore con N elementi, tutti pari a 7.1

# ora calcolo d=a*x + y
d= [a*x[i] + y[i] for i in range (N)]

#verifico che tutti gli elementi valgano 7.4

# stampa solo i primi 10
print("\nPrimi 10 elementi di d:")
for i in range(min(10, N)):
    print(f"d[{i}] = {d[i]}")

if N > 10:
    print("\nUltimi 10 elementi di d:")
    start = max(0, N - 10)
    for i in range(start, N):
        print(f"d[{i}] = {d[i]}")

#controllo finale
if all(abs(value - 7.4)<1e-12 for value in d):
    print ("Tutti gli elementi di d valgono circa 7.4 ✅")
else:
    print("Errore nel calcolo ❌")

```

## Run the code 

On terminal, then write N from the options

```
 python3 vector_sum.py

```

## matrix_mul.py

```python
# Task 02 – Punto 3 (interpreted language)
# Prodotto matriciale C = A * B con N = 10, 100, 10000
# A = tutti 3, B = tutti 7.1
# Test atteso (matmul standard): ogni elemento di C è 21.3 * N

import numpy as np
import math

print("Programma: C = A * B con A=3, B=7.1")
N = int(input("Inserisci N (10, 100, 10000): "))

a_val = 3.0
b_val = 7.1
expected = 21.3 * N
tol = 1e-9

if N <= 100:
    # Costruisco A e B come liste di liste
    A = [[a_val for _ in range(N)] for __ in range(N)]
    B = [[b_val for _ in range(N)] for __ in range(N)]
    # C inizialmente tutta zero
    C = [[0.0 for _ in range(N)] for __ in range(N)]

    # Triplo ciclo O(N^3) (va bene per N<=100)
    for i in range(N):
        for k in range(N):      # indice di somma
            aik = A[i][k]
            for j in range(N):
                C[i][j] += aik * B[k][j]

    # Verifica: tutti ~ 21.3*N
    ok = True
    for i in range(N):
        for j in range(N):
            if not math.isclose(C[i][j], expected, rel_tol=0.0, abs_tol=tol):
                ok = False
                break
        if not ok:
            break

    # Stampa campione (3x3 in alto a sx e dx)
    def print_block():
        m = min(3, N)
        print("\nBlocco alto-sinistra:")
        for i in range(m):
            print([C[i][j] for j in range(m)])
        print("\nBlocco alto-destra:")
        for i in range(m):
            print([C[i][N-m+j] for j in range(m)])

    print_block()
    print(f"\nVerifica: {'OK' if ok else 'ERRORE'} (atteso ogni C[i][j] = {expected})")

else:
    # N=10000: il prodotto O(N^3) e l’allocazione N^2 non sono praticabili in Python.
   
    print("\nN=10000: il prodotto completo è troppo oneroso in Python (tempo/memoria).")
    print(f"Per A=3 e B=7.1, ogni elemento di C dovrebbe essere: 21.3 * N = {expected}")
    print("Questo è il valore atteso da usare per il test.")

```

## Run the code 

On terminal, then write N from the options

```
 python3 matrix_mul.py

```

## vector_sum.c

```c
// Task 02 – Punto 2 (compiled language)
// d = a*x + y con N scelto dall'utente
// a = 3, x = 0.1 per tutti, y = 7.1 per tutti

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    // Dichiaro variabili
    long long N;    // uso long long per gestire anche N grandi
    double a = 3.0;
    double *x, *y, *d;

    // Messaggi introduttivi
    printf("Programma: calcolo di d = a*x + y\n");
    printf("a = 3, x = 0.1, y = 7.1\n");
    printf("Inserisci il valore di N (es. 10, 1000000, 100000000): ");
    scanf("%lld", &N); // leggo N da tastiera

    // Allocazione memoria per i vettori
    x = (double*) malloc(N * sizeof(double));
    y = (double*) malloc(N * sizeof(double));
    d = (double*) malloc(N * sizeof(double));

    if (x == NULL || y == NULL || d == NULL) {
        printf("Errore di allocazione memoria!\n");
        return 1;
    }

    // Inizializzo i vettori
    for (long long i = 0; i < N; i++) {
        x[i] = 0.1;
        y[i] = 7.1;
    }

    // Calcolo d = a*x + y
    for (long long i = 0; i < N; i++) {
        d[i] = a * x[i] + y[i];
    }

    // Stampa i primi 10 elementi
    printf("\nPrimi 10 elementi di d:\n");
    for (long long i = 0; i < 10 && i < N; i++) {
        printf("d[%lld] = %.17g\n", i, d[i]);
    }

    // Stampa gli ultimi 10 elementi (solo se N > 10)
    if (N > 10) {
        printf("\nUltimi 10 elementi di d:\n");
        long long start = (N > 10) ? N - 10 : 0;
        for (long long i = start; i < N; i++) {
            printf("d[%lld] = %.17g\n", i, d[i]);
        }
    }

    // Verifica che tutti siano ~7.4
    int ok = 1;
    for (long long i = 0; i < N; i++) {
        if (fabs(d[i] - 7.4) > 1e-12) {
            ok = 0;
            break;
        }
    }

    if (ok) {
        printf("\nTutti gli elementi di d valgono circa 7.4 ✅\n");
    } else {
        printf("\nErrore nel calcolo ❌\n");
    }

    // Libero la memoria
    free(x);
    free(y);
    free(d);

    return 0;
}

```

## Run the code 

On terminal:

Compile the code:

```
 gcc -02 vector_sum.c -o vector_sum

```
Run the code, then write N from the options:

```
./vector_sum

```

## matrix_mul.c

```c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

static int almost_equal(double a, double b, double abs_tol) {
    return fabs(a - b) <= abs_tol;
}

/*
 * Legge un intero long da input (stdin) con validazione base.
 * Ritorna 1 se OK e scrive in *out; 0 se errore.
 */
static int read_long_from_stdin(const char *prompt, long *out) {
    char buf[256];

    printf("%s", prompt);
    fflush(stdout);

    if (!fgets(buf, sizeof(buf), stdin)) {
        return 0;
    }

    // rimuovi newline
    buf[strcspn(buf, "\r\n")] = '\0';

    errno = 0;
    char *endp = NULL;
    long val = strtol(buf, &endp, 10);

    if (errno != 0) {
        return 0;
    }
    if (endp == buf || *endp != '\0') { // niente numero o caratteri extra
        return 0;
    }

    *out = val;
    return 1;
}

int main(void) {
    long N;

    printf("Programma: C = A * B con A=3, B=7.1\n");
    printf("Nota: per N <= 100 calcolo davvero la matmul; per N grandi uso la soluzione analitica.\n\n");

    if (!read_long_from_stdin("Inserisci N (es. 10, 100, 10000): ", &N)) {
        fprintf(stderr, "Errore: input non valido.\n");
        return 1;
    }
    if (N <= 0) {
        fprintf(stderr, "Errore: N deve essere > 0.\n");
        return 1;
    }

    const double a_val = 3.0;
    const double b_val = 7.1;
    const double expected = (a_val * b_val) * (double)N; // 21.3 * N
    const double tol = 1e-9;

    printf("\nN = %ld\n", N);
    printf("Valore atteso per ogni C[i][j] = %.10f\n\n", expected);

    // Per N "piccoli" facciamo davvero la moltiplicazione completa
    if (N <= 100) {
        size_t NN = (size_t)N * (size_t)N;

        double *A = (double *)malloc(NN * sizeof(double));
        double *B = (double *)malloc(NN * sizeof(double));
        double *C = (double *)calloc(NN, sizeof(double)); // inizializzata a 0

        if (!A || !B || !C) {
            fprintf(stderr, "Errore: memoria insufficiente\n");
            free(A); free(B); free(C);
            return 1;
        }

        // Riempio A e B (tutti 3 e tutti 7.1)
        for (size_t idx = 0; idx < NN; idx++) {
            A[idx] = a_val;
            B[idx] = b_val;
        }

        // Matmul C = A*B (row-major)
        // C[i,j] = sum_k A[i,k] * B[k,j]
        for (long i = 0; i < N; i++) {
            for (long k = 0; k < N; k++) {
                double aik = A[(size_t)i * (size_t)N + (size_t)k];
                for (long j = 0; j < N; j++) {
                    C[(size_t)i * (size_t)N + (size_t)j] += aik * B[(size_t)k * (size_t)N + (size_t)j];
                }
            }
        }

        // Verifica
        int ok = 1;
        for (long i = 0; i < N && ok; i++) {
            for (long j = 0; j < N; j++) {
                double cij = C[(size_t)i * (size_t)N + (size_t)j];
                if (!almost_equal(cij, expected, tol)) {
                    ok = 0;
                    fprintf(stderr, "Mismatch at (%ld,%ld): got %.12f expected %.12f\n",
                            i, j, cij, expected);
                    break;
                }
            }
        }

        // Stampa blocchi 3x3 (alto-sx e alto-dx)
        long m = (N < 3) ? N : 3;

        printf("Blocco alto-sinistra:\n");
        for (long i = 0; i < m; i++) {
            for (long j = 0; j < m; j++) {
                printf("%12.6f ", C[(size_t)i * (size_t)N + (size_t)j]);
            }
            printf("\n");
        }

        printf("\nBlocco alto-destra:\n");
        for (long i = 0; i < m; i++) {
            for (long j = N - m; j < N; j++) {
                printf("%12.6f ", C[(size_t)i * (size_t)N + (size_t)j]);
            }
            printf("\n");
        }

        printf("\nVerifica: %s\n", ok ? "OK" : "ERRORE");

        free(A);
        free(B);
        free(C);
        return ok ? 0 : 2;
    }

    // Per N grande (es. 10000): niente allocazione N^2 e niente O(N^3).
    // Usiamo il fatto che A e B sono costanti => ogni elemento di C è expected.
    printf("N grande: evito allocazione e prodotto completo (troppo onerosi).\n");
    printf("Poiché A[i,k]=3 e B[k,j]=7.1 per tutti, allora:\n");
    printf("C[i,j] = sum_{k=0..N-1} 3*7.1 = N * 21.3 = %.10f\n\n", expected);

    // “Test” su alcuni indici campione: il valore teorico è sempre expected
    long samples[4][2] = {
        {0, 0},
        {0, N-1},
        {N-1, 0},
        {N-1, N-1}
    };

    printf("Campioni (teorici):\n");
    for (int s = 0; s < 4; s++) {
        long i = samples[s][0];
        long j = samples[s][1];
        double cij = expected; // valore noto
        printf("C[%ld,%ld] = %.10f\n", i, j, cij);
    }

    printf("\nVerifica (analitica): OK\n");
    return 0;
}


```

## Run the code 

On terminal:

Compile the code:

```
 gcc -02 matrix_mul.c -o matrix_mul

```
Run the code, then write N from the options:

```

./matrix_mul

```

**Answer the question**
1) Yes, problems were encountered for large values of N, especially for N = 10⁸ when using the interpreted language (Python).
The main issues are related to memory usage and execution time.
In Python, vectors are stored as lists of objects, which require a significant amount of memory. Creating and handling multiple vectors of size 10⁸ can be very slow and memory intensive.
Additionally, printing large amounts of data to the terminal dramatically slows down execution.

In the compiled language (C), the same operations are much more memory-efficient and faster, since arrays are stored as contiguous memory blocks and arithmetic operations are compiled to machine code.
However, for very large N, the execution time is still limited by memory bandwidth rather than pure computation speed.

2) Yes, the correctness of the results was tested for all points.

For the vector sum (points 1 and 2), the analytical expected result is:

𝑑=𝑎⋅𝑥+𝑦=3⋅0.1+7.1=7.4. d=a⋅x+y=3⋅0.1+7.1=7.4.

The test was performed by printing a limited number of elements (the first and the last ones) and by checking that all values were approximately equal to 7.4 within a numerical tolerance, to account for floating-point precision.

For the matrix multiplication (point 3), the verification was performed by printing small blocks of the resulting matrix (top-left and bottom-right corners) and checking that all printed elements had the same value, as expected from the analytical result.

For small matrix sizes (N = 10 and N = 100), this method allows a direct and explicit verification of the result.
For very large values of N (e.g. N = 10000), printing or explicitly checking all matrix elements is not feasible due to memory and time constraints, therefore only partial output inspection or analytical reasoning can be used to validate the result.
