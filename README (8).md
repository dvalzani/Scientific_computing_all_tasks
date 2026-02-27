# Task09

## 1) Implementation of DAXPY code: daxpy_serial.c

```c
#include <stdio.h>      // printf, scanf
#include <stdlib.h>     // malloc, free
#include <sys/time.h>   // gettimeofday

// Funzione per misurare il tempo in secondi (wall clock)
double wall_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1e-6;
}

int main() {

    long N;

    // Richiesta input all'utente
    printf("Inserisci N: ");
    if (scanf("%ld", &N) != 1 || N <= 0) {
        printf("Errore: N deve essere un intero positivo\n");
        return 1;
    }

    // Alloco memoria per i vettori x, y, d
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));

    if (x == NULL || y == NULL || d == NULL) {
        printf("Errore allocazione memoria\n");
        return 1;
    }

    // Inizializzo i vettori
    for (long i = 0; i < N; i++) {
        x[i] = 1.0;
        y[i] = 2.0;
    }

    // Misuro il tempo del ciclo DAXPY
    double t0 = wall_time();

    for (long i = 0; i < N; i++) {
        d[i] = x[i] + y[i];
    }

    double t1 = wall_time();

    // Controllo veloce correttezza
    printf("Esempio: d[0] = %f (atteso 3.0)\n", d[0]);

    // Stampo il tempo impiegato
    printf("Serial time: %.6f sec (N = %ld)\n", t1 - t0, N);

    // Libero la memoria
    free(x);
    free(y);
    free(d);

    return 0;
}
 ```
## Compile

```c
gcc daxpy_serial.c -o daxpy_serial
```

## Run, then choose N

```c
./daxpy_serial
```


## 2) Implementation of DAXPY code usinh MPI directives: daxpy_mpi.c

```c
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/*
  Task 09 - DAXPY (semplificato): d = x + y
  Versione MPI con input interattivo + warning per N piccoli.

  Uso:
    - Interattivo:  mpirun -np 4 ./daxpy_mpi
    - Non-interattivo: mpirun -np 4 ./daxpy_mpi 4000000

  Note:
    - N deve essere multiplo del numero di processi (size), per semplicità.
    - Per N piccoli, MPI può risultare "lento" perché domina l'overhead di comunicazione:
      il codice stampa un WARNING ma NON blocca.
*/

// Soglia sotto cui il timing MPI è spesso dominato da overhead (non un errore)
#define N_WARN 1000000L

// Massimo "sicuro" per evitare che il rank 0 allochi troppa RAM nel container
#define N_MAX 20000000L

static void print_mem_estimate(long N, int size) {
    // rank 0: 3*N double (x,y,d)
    // ogni rank: 3*(N/size) double (x_local,y_local,d_local)
    double bytes_rank0 = 3.0 * (double)N * 8.0;
    double bytes_local = 3.0 * ((double)N / (double)size) * 8.0;

    printf("Stima memoria (circa):\n");
    printf("  Rank 0: %.1f MB (x,y,d globali)\n", bytes_rank0 / (1024.0 * 1024.0));
    printf("  Ogni rank: +%.1f MB (x_local,y_local,d_local)\n", bytes_local / (1024.0 * 1024.0));
    printf("  Nota: stima minima, ci sono overhead del runtime MPI.\n");
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    long N = 0;

    // ------------------------------------------------------------
    // 1) Scelta di N: da argv oppure interattivo (solo rank 0)
    // ------------------------------------------------------------
    if (argc >= 2) {
        // Modalità non-interattiva
        N = atol(argv[1]);

        if (N <= 0) {
            if (rank == 0) printf("Errore: N deve essere > 0\n");
            MPI_Finalize();
            return 1;
        }

        if (N > N_MAX) {
            if (rank == 0) {
                printf("Errore: N=%ld supera il massimo di sicurezza N_MAX=%ld\n", N, N_MAX);
                printf("Riduci N per evitare problemi di memoria nel container.\n");
            }
            MPI_Finalize();
            return 1;
        }

        if (N % size != 0) {
            if (rank == 0) {
                printf("Errore: N (%ld) deve essere multiplo del numero di processi (%d)\n", N, size);
            }
            MPI_Finalize();
            return 1;
        }

    } else {
        // Modalità interattiva: solo rank 0 chiede N, poi broadcast
        if (rank == 0) {
            while (1) {
                printf("Inserisci N (1 .. %ld) multiplo di %d: ", N_MAX, size);
                fflush(stdout);

                if (scanf("%ld", &N) != 1) {
                    // input non valido: pulisco stdin e riprovo
                    int c;
                    while ((c = getchar()) != '\n' && c != EOF) {}
                    printf("Input non valido. Riprova.\n");
                    continue;
                }

                if (N <= 0) {
                    printf("N deve essere > 0. Riprova.\n");
                    continue;
                }

                if (N > N_MAX) {
                    printf("N troppo grande (max %ld). Riprova.\n", N_MAX);
                    continue;
                }

                if (N % size != 0) {
                    printf("N deve essere multiplo di %d. Riprova.\n", size);
                    continue;
                }

                // (opzionale ma utile) stima memoria
                print_mem_estimate(N, size);

                break;
            }
        }

        MPI_Bcast(&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    }

    // ------------------------------------------------------------
    // 2) WARNING per N piccoli (non blocca)
    // ------------------------------------------------------------
    if (rank == 0 && N < N_WARN) {
        printf("\n[WARNING]\n");
        printf("  N = %ld e' relativamente piccolo per benchmark MPI.\n", N);
        printf("  Il tempo puo' essere dominato da MPI_Init / comunicazione / sincronizzazioni,\n");
        printf("  piu' che dal kernel DAXPY. Questo e' normale e NON indica un errore.\n\n");
    }

    long local_N = N / size;

    // ------------------------------------------------------------
    // 3) Allocazioni globali (solo rank 0)
    // ------------------------------------------------------------
    double *x = NULL;
    double *y = NULL;
    double *d = NULL;

    if (rank == 0) {
        x = (double*) malloc(N * sizeof(double));
        y = (double*) malloc(N * sizeof(double));
        d = (double*) malloc(N * sizeof(double));

        if (x == NULL || y == NULL || d == NULL) {
            printf("Errore: allocazione memoria globale fallita nel rank 0 (N=%ld)\n", N);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (long i = 0; i < N; i++) {
            x[i] = 1.0;
            y[i] = 2.0;
        }
    }

    // ------------------------------------------------------------
    // 4) Allocazioni locali (tutti i rank)
    // ------------------------------------------------------------
    double *x_local = (double*) malloc(local_N * sizeof(double));
    double *y_local = (double*) malloc(local_N * sizeof(double));
    double *d_local = (double*) malloc(local_N * sizeof(double));

    if (x_local == NULL || y_local == NULL || d_local == NULL) {
        printf("Rank %d: errore allocazione memoria locale (local_N=%ld)\n", rank, local_N);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // ------------------------------------------------------------
    // 5) Scatter dati
    // ------------------------------------------------------------
    MPI_Scatter(x, local_N, MPI_DOUBLE,
                x_local, local_N, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    MPI_Scatter(y, local_N, MPI_DOUBLE,
                y_local, local_N, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // ------------------------------------------------------------
    // 6) Timing (kernel only): misura SOLO il loop DAXPY locale
    // ------------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (long i = 0; i < local_N; i++) {
        d_local[i] = x_local[i] + y_local[i];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    // ------------------------------------------------------------
    // 7) Gather risultati nel rank 0
    // ------------------------------------------------------------
    MPI_Gather(d_local, local_N, MPI_DOUBLE,
               d,       local_N, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // ------------------------------------------------------------
    // 8) Reduction per controllo correttezza (somma globale)
    // ------------------------------------------------------------
    double local_sum = 0.0;
    for (long i = 0; i < local_N; i++) {
        local_sum += d_local[i];
    }

    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // ------------------------------------------------------------
    // 9) Output finale (solo rank 0)
    // ------------------------------------------------------------
    if (rank == 0) {
        printf("Esempio: d[0] = %f (atteso 3.0)\n", d[0]);
        printf("MPI time (kernel only): %.6f sec (N = %ld, processes = %d)\n", t1 - t0, N, size);
        printf("Somma globale di d = %.10f (atteso ~ %.10f)\n", global_sum, 3.0 * (double)N);

        free(x);
        free(y);
        free(d);
    }

    free(x_local);
    free(y_local);
    free(d_local);

    MPI_Finalize();
    return 0;
}

```
## Compile (Makefile at the end of this repositories)

```c
make
```
## Run, then choose N

```c
make run-mpi
```

## 3) Optional

```c

#include <stdio.h>      // printf, scanf
#include <stdlib.h>     // malloc, free
#include <sys/time.h>   // gettimeofday
#include <math.h>       // fabs
#include <omp.h>        // OpenMP

// Wall-clock time in seconds
static double wall_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1e-6;
}

int main() {
    long N;
    int threads;

    // ---- Input interattivo: N ----
    printf("Inserisci N: ");
    fflush(stdout);
    if (scanf("%ld", &N) != 1 || N <= 0) {
        printf("Errore: N deve essere un intero positivo\n");
        return 1;
    }

    // ---- Input interattivo: threads ----
    printf("Inserisci numero di threads OpenMP: ");
    fflush(stdout);
    if (scanf("%d", &threads) != 1 || threads <= 0) {
        printf("Errore: threads deve essere un intero positivo\n");
        return 1;
    }

    // Warning: per N piccoli i tempi OpenMP possono essere dominati dall'overhead
    if (N < 1000000L) {
        printf("\n[WARNING]\n");
        printf("  N = %ld e' relativamente piccolo per benchmark OpenMP.\n", N);
        printf("  Il tempo puo' essere dominato da overhead di threads/scheduling.\n");
        printf("  Questo e' normale e NON indica un errore.\n\n");
    }

    // ---- Allocazioni ----
    double *x = (double*) malloc((size_t)N * sizeof(double));
    double *y = (double*) malloc((size_t)N * sizeof(double));
    double *d = (double*) malloc((size_t)N * sizeof(double));

    if (x == NULL || y == NULL || d == NULL) {
        printf("Errore allocazione memoria (N=%ld)\n", N);
        free(x); free(y); free(d);
        return 1;
    }

    // ---- Inizializzazione ----
    for (long i = 0; i < N; i++) {
        x[i] = 1.0;
        y[i] = 2.0;
    }

    // Imposto numero di threads
    omp_set_num_threads(threads);

    // ============================================================
    // 1) DAXPY OpenMP (punto 1): d[i] = x[i] + y[i]
    // ============================================================
    double t0 = wall_time();

    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        d[i] = x[i] + y[i];
    }

    double t1 = wall_time();

    // ============================================================
    // 3) (Optional) Reduction: total_sum = sum(d)
    //    - seriale (baseline)
    //    - OpenMP reduction (quello richiesto dalla task)
    // ============================================================

    // Somma seriale (baseline)
    double sum_serial = 0.0;
    for (long i = 0; i < N; i++) {
        sum_serial += d[i];
    }

    // Somma OpenMP con reduction (punto 3!)
    double sum_omp = 0.0;
    double t2 = wall_time();

    #pragma omp parallel for reduction(+:sum_omp)
    for (long i = 0; i < N; i++) {
        sum_omp += d[i];
    }

    double t3 = wall_time();

    // Atteso: d[i] = 3.0, quindi sum attesa = 3.0 * N
    double expected = 3.0 * (double)N;

    // Errori relativi (evitiamo divisione per zero)
    double rel_err_serial = fabs(sum_serial - expected) / (fabs(expected) + 1e-30);
    double rel_err_omp    = fabs(sum_omp    - expected) / (fabs(expected) + 1e-30);
    double rel_diff       = fabs(sum_omp - sum_serial) / (fabs(sum_serial) + 1e-30);

    // ---- Output ----
    printf("Esempio: d[0] = %f (atteso 3.0)\n", d[0]);
    printf("OpenMP time (DAXPY): %.6f sec (N = %ld, threads = %d)\n", t1 - t0, N, threads);

    printf("\n--- Reduction check (Task 09 - punto 3 opzionale) ---\n");
    printf("Somma attesa      = %.10f\n", expected);
    printf("Somma seriale     = %.10f (rel_err = %.3e)\n", sum_serial, rel_err_serial);
    printf("Somma OpenMP red. = %.10f (rel_err = %.3e)\n", sum_omp, rel_err_omp);
    printf("Diff OpenMP-serial (rel) = %.3e\n", rel_diff);
    printf("OpenMP time (reduction): %.6f sec\n", t3 - t2);

    // ---- Free ----
    free(x);
    free(y);
    free(d);

    return 0;
}

```

## Compile (Makefile at the end of this repositories)

```c
make
```
## Run, then choose N and number of threads (no more or equal than 4)

```c
./daxpy_openmp
```

## 3) Makefile

```c

CC=gcc
CFLAGS=-O3

MPICC=/usr/lib64/openmpi/bin/mpicc
MPIRUN=/usr/lib64/openmpi/bin/mpirun

# OpenMPI settings that avoid container issues (OFI/UCX)
MPI_ENV=OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        OMPI_MCA_pml=ob1 OMPI_MCA_osc=pt2pt OMPI_MCA_mtl=^ofi OMPI_MCA_btl=self,vader \
        LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$$LD_LIBRARY_PATH

all: daxpy_serial daxpy_openmp daxpy_mpi

daxpy_serial: daxpy_serial.c
	$(CC) $(CFLAGS) $< -o $@

daxpy_openmp: daxpy_openmp.c
	$(CC) $(CFLAGS) -fopenmp $< -o $@

daxpy_mpi: daxpy_mpi.c
	$(MPI_ENV) $(MPICC) $(CFLAGS) $< -o $@

run-serial: daxpy_serial
	./daxpy_serial

run-omp: daxpy_openmp
	./daxpy_openmp

run-mpi: daxpy_mpi
	$(MPI_ENV) $(MPIRUN) -np 4 ./daxpy_mpi

clean:
	rm -f daxpy_serial daxpy_openmp daxpy_mpi

```
