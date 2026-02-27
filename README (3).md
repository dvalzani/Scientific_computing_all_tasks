# Task 03 – Basis of programming

## 1 `generate_vectors.c`
```c
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int N;
    char prefix[128];

    // Caso 1: argomenti da linea di comando
    if (argc >= 3) {
        N = atoi(argv[1]);
        snprintf(prefix, sizeof(prefix), "%s", argv[2]);
    } 
    // Caso 2: input interattivo
    else {
        printf("Enter N (vector size): ");
        if (scanf("%d", &N) != 1) {
            printf("Error: invalid N.\n");
            return 1;
        }

        printf("Enter output prefix (e.g., ./vector_): ");
        if (scanf("%127s", prefix) != 1) {
            printf("Error: invalid prefix.\n");
            return 1;
        }
    }

    if (N <= 0) {
        printf("Error: N must be > 0.\n");
        return 1;
    }

    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));

    if (x == NULL || y == NULL) {
        perror("Error allocating memory");
        free(x);
        free(y);
        return 1;
    }

    // Fill vectors
    for (int i = 0; i < N; i++) {
        x[i] = 0.1;
        y[i] = 7.1;
    }

    // Output file names
    char filename_x[256], filename_y[256];
    snprintf(filename_x, sizeof(filename_x), "%sN%d_x.dat", prefix, N);
    snprintf(filename_y, sizeof(filename_y), "%sN%d_y.dat", prefix, N);

    FILE *fx = fopen(filename_x, "w");
    FILE *fy = fopen(filename_y, "w");

    if (fx == NULL || fy == NULL) {
        perror("Error opening output files");
        if (fx) fclose(fx);
        if (fy) fclose(fy);
        free(x);
        free(y);
        return 1;
    }

    // Write files
    for (int i = 0; i < N; i++) {
        fprintf(fx, "%.2f\n", x[i]);
        fprintf(fy, "%.2f\n", y[i]);
    }

    fclose(fx);
    fclose(fy);
    free(x);
    free(y);

    printf("Files generated:\n  %s\n  %s\n", filename_x, filename_y);
    return 0;
}

```
### 1.1 Compilation and Execution

**To compile (on terminal)**:

```

gcc generate_vectors.c -o generate_vectors

```

**To run the program, then choose N and select output name**
```
./generate_vectors 

```

**After execution, the following files will be generated (eg "vector" name and 10 for N is chosen):**

```
vectorN10_x.dat
vectorN10_y.dat

```


## 2 `compute_ax_plus_y.c`

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static void trim_in_place(char *s) {
    // trim left
    char *p = s;
    while (*p && isspace((unsigned char)*p)) p++;
    if (p != s) memmove(s, p, strlen(p) + 1);

    // trim right
    size_t len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) {
        s[len - 1] = '\0';
        len--;
    }
}

static int file_exists(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static int read_kv(FILE *f, const char *key, char *out, size_t out_sz) {
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        trim_in_place(line);
        if (line[0] == '\0' || line[0] == '#') continue;

        char *eq = strchr(line, '=');
        if (!eq) continue;

        *eq = '\0';
        char left[256];
        snprintf(left, sizeof(left), "%s", line);
        trim_in_place(left);

        char right[256];
        snprintf(right, sizeof(right), "%s", eq + 1);
        trim_in_place(right);

        if (strcmp(left, key) == 0) {
            snprintf(out, out_sz, "%s", right);
            return 1;
        }
    }
    return 0;
}

static int write_config(const char *path,
                        const char *file_x,
                        const char *file_y,
                        int N,
                        double a,
                        const char *prefix_output) {
    FILE *fc = fopen(path, "w");
    if (!fc) return 0;

    fprintf(fc, "file_x=%s\n", file_x);
    fprintf(fc, "file_y=%s\n", file_y);
    fprintf(fc, "N=%d\n", N);
    fprintf(fc, "a=%.17g\n", a);
    fprintf(fc, "prefix_output=%s\n", prefix_output);

    fclose(fc);
    return 1;
}

int main(void) {
    char file_x[256] = "";
    char file_y[256] = "";
    char prefix_output[256] = "vector_";
    int N = 0;
    double a = 0.0;

    const char *config_name = "config.txt";

    // 1) If config.txt does not exist, ask user and create it
    if (!file_exists(config_name)) {
        printf("config.txt not found. Let's create it.\n");

        printf("Enter filename for x (e.g., vector_N10_x.dat): ");
        if (scanf("%255s", file_x) != 1) return 1;

        printf("Enter filename for y (e.g., vector_N10_y.dat): ");
        if (scanf("%255s", file_y) != 1) return 1;

        printf("Enter N (vector size): ");
        if (scanf("%d", &N) != 1) return 1;

        printf("Enter a (scalar): ");
        if (scanf("%lf", &a) != 1) return 1;

        printf("Enter output prefix (e.g., vector_): ");
        if (scanf("%255s", prefix_output) != 1) return 1;

        if (N <= 0) {
            printf("Error: N must be > 0\n");
            return 1;
        }

        if (!write_config(config_name, file_x, file_y, N, a, prefix_output)) {
            printf("Error: cannot write %s\n", config_name);
            return 1;
        }

        printf("Created %s\n\n", config_name);
    }

    // 2) Read config.txt
    FILE *fconf = fopen(config_name, "r");
    if (!fconf) {
        printf("Error: cannot open %s\n", config_name);
        return 1;
    }

    // We read each key by rewinding file each time (simple & clear).
    char tmp[256];

    rewind(fconf);
    if (!read_kv(fconf, "file_x", tmp, sizeof(tmp))) { printf("Missing file_x\n"); fclose(fconf); return 1; }
    snprintf(file_x, sizeof(file_x), "%s", tmp);

    rewind(fconf);
    if (!read_kv(fconf, "file_y", tmp, sizeof(tmp))) { printf("Missing file_y\n"); fclose(fconf); return 1; }
    snprintf(file_y, sizeof(file_y), "%s", tmp);

    rewind(fconf);
    if (!read_kv(fconf, "N", tmp, sizeof(tmp))) { printf("Missing N\n"); fclose(fconf); return 1; }
    N = atoi(tmp);

    rewind(fconf);
    if (!read_kv(fconf, "a", tmp, sizeof(tmp))) { printf("Missing a\n"); fclose(fconf); return 1; }
    a = atof(tmp);

    rewind(fconf);
    if (!read_kv(fconf, "prefix_output", tmp, sizeof(tmp))) { printf("Missing prefix_output\n"); fclose(fconf); return 1; }
    snprintf(prefix_output, sizeof(prefix_output), "%s", tmp);

    fclose(fconf);

    if (N <= 0) {
        printf("Error: invalid N (%d)\n", N);
        return 1;
    }

    printf("Using config:\n");
    printf("  file_x=%s\n", file_x);
    printf("  file_y=%s\n", file_y);
    printf("  N=%d\n", N);
    printf("  a=%.6f\n", a);
    printf("  prefix_output=%s\n\n", prefix_output);

    // 3) Allocate arrays
    double *x = malloc((size_t)N * sizeof(double));
    double *y = malloc((size_t)N * sizeof(double));
    double *d = malloc((size_t)N * sizeof(double));
    if (!x || !y || !d) {
        printf("Error: memory allocation failed\n");
        free(x); free(y); free(d);
        return 1;
    }

    // 4) Read input files
    FILE *fx = fopen(file_x, "r");
    FILE *fy = fopen(file_y, "r");
    if (!fx || !fy) {
        printf("Error: cannot open input files\n");
        if (fx) fclose(fx);
        if (fy) fclose(fy);
        free(x); free(y); free(d);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        if (fscanf(fx, "%lf", &x[i]) != 1) {
            printf("Error reading x at i=%d\n", i);
            fclose(fx); fclose(fy);
            free(x); free(y); free(d);
            return 1;
        }
        if (fscanf(fy, "%lf", &y[i]) != 1) {
            printf("Error reading y at i=%d\n", i);
            fclose(fx); fclose(fy);
            free(x); free(y); free(d);
            return 1;
        }
    }

    fclose(fx);
    fclose(fy);

    // 5) Compute d = a*x + y
    for (int i = 0; i < N; i++) {
        d[i] = a * x[i] + y[i];
    }

    // 6) Write output
    char filename_d[512];
    snprintf(filename_d, sizeof(filename_d), "%sN%d_d.dat", prefix_output, N);

    FILE *fd = fopen(filename_d, "w");
    if (!fd) {
        printf("Error: cannot open output file\n");
        free(x); free(y); free(d);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fprintf(fd, "%.2f\n", d[i]);
    }
    fclose(fd);

    free(x);
    free(y);
    free(d);

    printf("Done. Result saved in %s\n", filename_d);
    return 0;
}

```
### 2.1 Compilation and Execution

**To compile (on terminal)**:

```
gcc compute_ax_plus_y.c -o compute_ax_plus_y
```

**To run the program**
```
./compute_ax_plus_y
```

**Then, to create the config.txt and run correctly the script, write the command as follow (note that I used the previous vector files, you can choose other vector set with different N and name)**
```
config.txt not found. Let's create it.
Enter filename for x (e.g., vector_N10_x.dat): vectorN10_x.dat
Enter filename for y (e.g., vector_N10_y.dat): vectorN10_y.dat
Enter N (vector size): 10
Enter a (scalar): 3.00
Enter output prefix (e.g., vector_): daxpy
Created config.txt

Using config:
  file_x=vectorN10_x.dat
  file_y=vectorN10_y.dat
  N=10
  a=3.000000
  prefix_output=daxpy

Done. Result saved in daxpyN10_d.dat

```

**If you want to do other calculations, be sure to delete the last config.txt before run again compute_ax_plus_y.c**

## 3 Makefile

To compile both programs (`generate_vectors.c` and `compute_ax_plus_y.c`) we can use the following `Makefile`:

```makefile
CC     = gcc
CFLAGS = -Wall -Wextra -O2

# Main target: compile both programs
all: generate_vectors compute_ax_plus_y

# Compile generate_vectors
generate_vectors: generate_vectors.c
	$(CC) $(CFLAGS) generate_vectors.c -o generate_vectors

# Compile compute_ax_plus_y
compute_ax_plus_y: compute_ax_plus_y.c
	$(CC) $(CFLAGS) compute_ax_plus_y.c -o compute_ax_plus_y

# Clean all generated files
clean:
	rm -f generate_vectors compute_ax_plus_y *.o vector_*.dat
```
**To compile all the scripts, on terminal**

```c
make
```

## 4 Code for HDF5 files:

Before modify our code to use and read hdf5 files, we have to update the Makefile for the points 4 and 5 of the task03:

```makefile
CC       = gcc
CFLAGS   = -Wall -Wextra -O2
GSL_LIBS = -lgsl -lgslcblas -lm

# ------------------------------------------------------------
# HDF5: no pkg-config available in this environment.
# Use detected include path + prefer libhdf5_openmpi if present.
# ------------------------------------------------------------
HDF5_INC  = -I/usr/include/openmpi-x86_64
HDF5_LIBS = $(shell ldconfig -p 2>/dev/null | grep -q libhdf5_openmpi && echo -lhdf5_openmpi || echo -lhdf5)


# ------------------------------------------------------------
# Targets
# ------------------------------------------------------------
all: generate_vectors compute_ax_plus_y compute_ax_plus_y_gsl \
     generate_vectors_hdf5 compute_ax_plus_y_hdf5

generate_vectors: generate_vectors.c
	$(CC) $(CFLAGS) $< -o $@

compute_ax_plus_y: compute_ax_plus_y.c
	$(CC) $(CFLAGS) $< -o $@

compute_ax_plus_y_gsl: compute_ax_plus_y_gsl.c
	$(CC) $(CFLAGS) $< -o $@ $(GSL_LIBS)

generate_vectors_hdf5: generate_vectors_hdf5.c
	$(CC) $(CFLAGS) $(HDF5_INC) $< -o $@ $(HDF5_LIBS)

compute_ax_plus_y_hdf5: compute_ax_plus_y_hdf5.c
	$(CC) $(CFLAGS) $(HDF5_INC) $< -o $@ $(HDF5_LIBS)


clean:
	rm -f generate_vectors compute_ax_plus_y compute_ax_plus_y_gsl \
	      generate_vectors_hdf5 compute_ax_plus_y_hdf5 \
	      *.o vector_*.dat *.h5

```

### 4.1 generate vectors in hdf5 file format

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"

int main(int argc, char *argv[]) {

    int N = 0;
    char filename[256];

    // --- Mode selection: CLI or interactive ---
    if (argc >= 3) {
        // CLI mode: ./generate_vectors_hdf5 N output.h5
        N = atoi(argv[1]);
        snprintf(filename, sizeof(filename), "%s", argv[2]);
    } else {
        // Interactive mode: ./generate_vectors_hdf5
        printf("Enter N (vector size): ");
        if (scanf("%d", &N) != 1) {
            printf("Error: invalid input for N.\n");
            return 1;
        }
        snprintf(filename, sizeof(filename), "vectors_N%d.h5", N);
    }

    if (N <= 0) {
        printf("Error: N must be > 0\n");
        return 1;
    }

    // Allocate vectors
    double *x = malloc((size_t)N * sizeof(double));
    double *y = malloc((size_t)N * sizeof(double));
    if (x == NULL || y == NULL) {
        perror("Error allocating memory");
        free(x);
        free(y);
        return 1;
    }

    // Fill vectors
    for (int i = 0; i < N; i++) {
        x[i] = 0.1;
        y[i] = 7.1;
    }

    // --- HDF5 writing ---
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        printf("Error creating HDF5 file %s\n", filename);
        free(x);
        free(y);
        return 1;
    }

    hsize_t dims[1] = { (hsize_t)N };
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    if (dataspace_id < 0) {
        printf("Error creating dataspace\n");
        H5Fclose(file_id);
        free(x);
        free(y);
        return 1;
    }

    // Dataset /x
    hid_t dset_x_id = H5Dcreate(file_id, "/x", H5T_NATIVE_DOUBLE,
                                dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset_x_id < 0) {
        printf("Error creating dataset /x\n");
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        free(x);
        free(y);
        return 1;
    }
    if (H5Dwrite(dset_x_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x) < 0) {
        printf("Error writing dataset /x\n");
        H5Dclose(dset_x_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        free(x);
        free(y);
        return 1;
    }
    H5Dclose(dset_x_id);

    // Dataset /y
    hid_t dset_y_id = H5Dcreate(file_id, "/y", H5T_NATIVE_DOUBLE,
                                dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset_y_id < 0) {
        printf("Error creating dataset /y\n");
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        free(x);
        free(y);
        return 1;
    }
    if (H5Dwrite(dset_y_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y) < 0) {
        printf("Error writing dataset /y\n");
        H5Dclose(dset_y_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        free(x);
        free(y);
        return 1;
    }
    H5Dclose(dset_y_id);

    H5Sclose(dataspace_id);
    H5Fclose(file_id);

    free(x);
    free(y);

    printf("HDF5 file generated: %s (datasets: /x, /y, size N = %d)\n", filename, N);
    return 0;
}

```

**On terminal, now choose the N value (e.g. N=50), a file .hdf5 is created as follow:**

```

./generate_vectors_hdf5 
Enter N (vector size): 50
HDF5 file generated: vectors_N50.h5 (datasets: /x, /y, size N = 50)

```

**Before perform d=ax+y, we need a config file. We can write it as follow and name it config_hdf5.txt:**

```config_hdf5.txt
file_h5=vectors_N50.h5
dataset_x=/x
dataset_y=/y
N=50
a=3.0
prefix_output=vector_
```

### 4.2 read info from HDF5 file and perform d=ax+y

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"

// Funzione per estrarre il valore dopo '=' e togliere newline
void leggi_valore(char *riga, char *valore) {
    char *uguale = strchr(riga, '=');
    if (uguale != NULL) {
        strcpy(valore, uguale + 1);
        // Rimuove newline
        valore[strcspn(valore, "\r\n")] = '\0';
        // Rimuove spazi iniziali
        while (*valore == ' ') valore++;
    }
}

int main(void) {
    char file_h5[128], dataset_x[128], dataset_y[128], prefix_output[128];
    int N = 0;
    double a = 0.0;

    FILE *fconf = fopen("config_hdf5.txt", "r");
    if (fconf == NULL) {
        printf("Error: cannot open config_hdf5.txt\n");
        return 1;
    }

    char riga[256], valore[128];

    // Parsing del file di configurazione
    while (fgets(riga, sizeof(riga), fconf)) {
        if (strstr(riga, "file_h5")) {
            leggi_valore(riga, file_h5);
        } else if (strstr(riga, "dataset_x")) {
            leggi_valore(riga, dataset_x);
        } else if (strstr(riga, "dataset_y")) {
            leggi_valore(riga, dataset_y);
        } else if (strstr(riga, "N")) {
            leggi_valore(riga, valore);
            N = atoi(valore);
        } else if (strstr(riga, "a")) {
            leggi_valore(riga, valore);
            a = atof(valore);
        } else if (strstr(riga, "prefix_output")) {
            leggi_valore(riga, prefix_output);
        }
    }

    fclose(fconf);

    if (N <= 0) {
        printf("Error: invalid N (%d)\n", N);
        return 1;
    }

    printf("HDF5 input file : %s\n", file_h5);
    printf("Dataset x       : %s\n", dataset_x);
    printf("Dataset y       : %s\n", dataset_y);
    printf("N               : %d\n", N);
    printf("a               : %.2f\n", a);
    printf("Output prefix   : %s\n", prefix_output);

    // Alloca i vettori
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));

    if (x == NULL || y == NULL || d == NULL) {
        printf("Error allocating memory.\n");
        return 1;
    }

    // --- Lettura da HDF5 ---

    // Apre il file HDF5 in sola lettura
    hid_t file_id = H5Fopen(file_h5, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        printf("Error opening HDF5 file %s\n", file_h5);
        free(x); free(y); free(d);
        return 1;
    }

    // Apre il dataset x
    hid_t dset_x_id = H5Dopen(file_id, dataset_x, H5P_DEFAULT);
    if (dset_x_id < 0) {
        printf("Error opening dataset %s\n", dataset_x);
        H5Fclose(file_id);
        free(x); free(y); free(d);
        return 1;
    }

    // Legge x
    if (H5Dread(dset_x_id, H5T_NATIVE_DOUBLE,
                H5S_ALL, H5S_ALL, H5P_DEFAULT, x) < 0) {
        printf("Error reading dataset %s\n", dataset_x);
        H5Dclose(dset_x_id);
        H5Fclose(file_id);
        free(x); free(y); free(d);
        return 1;
    }
    H5Dclose(dset_x_id);

    // Apre il dataset y
    hid_t dset_y_id = H5Dopen(file_id, dataset_y, H5P_DEFAULT);
    if (dset_y_id < 0) {
        printf("Error opening dataset %s\n", dataset_y);
        H5Fclose(file_id);
        free(x); free(y); free(d);
        return 1;
    }

    // Legge y
    if (H5Dread(dset_y_id, H5T_NATIVE_DOUBLE,
                H5S_ALL, H5S_ALL, H5P_DEFAULT, y) < 0) {
        printf("Error reading dataset %s\n", dataset_y);
        H5Dclose(dset_y_id);
        H5Fclose(file_id);
        free(x); free(y); free(d);
        return 1;
    }
    H5Dclose(dset_y_id);

    H5Fclose(file_id);

    // --- Calcolo d = a*x + y ---
    for (int i = 0; i < N; i++) {
        d[i] = a * x[i] + y[i];
    }

    // --- Scrittura del risultato in HDF5 ---
    char filename_d[256];
    snprintf(filename_d, sizeof(filename_d), "%sN%d_d_hdf5.h5", prefix_output, N);

    hid_t file_out_id = H5Fcreate(filename_d, H5F_ACC_TRUNC,
                                  H5P_DEFAULT, H5P_DEFAULT);
    if (file_out_id < 0) {
        printf("Error creating output HDF5 file %s\n", filename_d);
        free(x); free(y); free(d);
        return 1;
    }

    hsize_t dims[1] = { (hsize_t) N };
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

    if (dataspace_id < 0) {
        printf("Error creating dataspace for output\n");
        H5Fclose(file_out_id);
        free(x); free(y); free(d);
        return 1;
    }

    hid_t dset_d_id = H5Dcreate(file_out_id, "/d",
                                H5T_NATIVE_DOUBLE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset_d_id < 0) {
        printf("Error creating dataset /d\n");
        H5Sclose(dataspace_id);
        H5Fclose(file_out_id);
        free(x); free(y); free(d);
        return 1;
    }

    if (H5Dwrite(dset_d_id, H5T_NATIVE_DOUBLE,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT, d) < 0) {
        printf("Error writing dataset /d\n");
        H5Dclose(dset_d_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_out_id);
        free(x); free(y); free(d);
        return 1;
    }

    H5Dclose(dset_d_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_out_id);

    printf("Operation completed (HDF5). Result saved in %s (dataset /d)\n",
           filename_d);

    free(x);
    free(y);
    free(d);

    return 0;
}


```

**Now, to calculate d=ax+y:**

```
./compute_ax_plus_y_hdf5

```

**and to see the result:**
```
h5dump -d /d vector_N50_d_hdf5.h5

```

## 5 Code for gsl_vector_axpby (you need the configfile.txt of the point1):

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>

// Funzione per estrarre il valore dopo '='
void leggi_valore(char *riga, char *valore) {
    char *uguale = strchr(riga, '=');
    if (uguale != NULL) {
        strcpy(valore, uguale + 1);
        valore[strcspn(valore, "\r\n")] = '\0';  // rimuove \n o \r\n
        while (*valore == ' ') valore++;         // rimuove spazi iniziali
    }
}

int main(void) {
    char file_x[128], file_y[128], prefix_output[128];
    int N = 0;
    double a = 0.0;

    FILE *fconf = fopen("config.txt", "r");
    if (fconf == NULL) {
        printf("Errore: impossibile aprire config.txt\n");
        return 1;
    }

    char riga[256], valore[128];

    // Parsing del file di configurazione
    while (fgets(riga, sizeof(riga), fconf)) {
        if (strstr(riga, "file_x")) {
            leggi_valore(riga, file_x);
        } else if (strstr(riga, "file_y")) {
            leggi_valore(riga, file_y);
        } else if (strstr(riga, "N")) {
            leggi_valore(riga, valore);
            N = atoi(valore);
        } else if (strstr(riga, "a")) {
            leggi_valore(riga, valore);
            a = atof(valore);
        } else if (strstr(riga, "prefix_output")) {
            leggi_valore(riga, prefix_output);
        }
    }

    fclose(fconf);

    if (N <= 0) {
        printf("Errore: N non valido (%d)\n", N);
        return 1;
    }

    printf("File x: %s\n", file_x);
    printf("File y: %s\n", file_y);
    printf("N: %d\n", N);
    printf("a: %.2f\n", a);
    printf("Output prefix: %s\n", prefix_output);

    // Alloca vettori GSL
    gsl_vector *vx = gsl_vector_alloc(N);
    gsl_vector *vy = gsl_vector_alloc(N);

    if (vx == NULL || vy == NULL) {
        printf("Errore allocazione vettori GSL.\n");
        return 1;
    }

    // Apre i file di input
    FILE *fx = fopen(file_x, "r");
    FILE *fy = fopen(file_y, "r");

    if (fx == NULL || fy == NULL) {
        printf("Errore apertura file x o y.\n");
        gsl_vector_free(vx);
        gsl_vector_free(vy);
        return 1;
    }

    // Legge i dati e li inserisce nei vettori GSL
    for (int i = 0; i < N; i++) {
        double valx, valy;
        if (fscanf(fx, "%lf", &valx) != 1) {
            printf("Errore lettura x[%d]\n", i);
            fclose(fx);
            fclose(fy);
            gsl_vector_free(vx);
            gsl_vector_free(vy);
            return 1;
        }
        if (fscanf(fy, "%lf", &valy) != 1) {
            printf("Errore lettura y[%d]\n", i);
            fclose(fx);
            fclose(fy);
            gsl_vector_free(vx);
            gsl_vector_free(vy);
            return 1;
        }
        gsl_vector_set(vx, i, valx);
        gsl_vector_set(vy, i, valy);
    }

    fclose(fx);
    fclose(fy);

    // Calcolo con GSL: vy = a*vx + 1.0*vy  → vy contiene d = a*x + y
    gsl_vector_axpby(a, vx, 1.0, vy);

    // Costruzione del nome del file di output
    char filename_d[256];
    snprintf(filename_d, sizeof(filename_d), "%sN%d_d_gsl.dat", prefix_output, N);

    FILE *fd = fopen(filename_d, "w");
    if (fd == NULL) {
        printf("Errore apertura file di output.\n");
        gsl_vector_free(vx);
        gsl_vector_free(vy);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        double di = gsl_vector_get(vy, i);
        fprintf(fd, "%.2f\n", di);
    }

    fclose(fd);

    gsl_vector_free(vx);
    gsl_vector_free(vy);

    printf("Operazione completata (GSL). Risultato salvato in %s\n", filename_d);
    return 0;
}
```
**with the config.txt of the point 1, now on terminal**

```c
./compute_ax_plus_y_gsl
```

**The output**
```
File x: vectorN10_x.dat
File y: vectorN10_y.dat
N: 10
a: 0.00
Output prefix: 
Operazione completata (GSL). Risultato salvato in N10_d_gsl.dat
```
