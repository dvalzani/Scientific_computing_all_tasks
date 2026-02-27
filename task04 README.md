## Given the folloving function

$f(x) = e^x \cos(x)$

We write a program in compiled language (c) that:

### INPUT

- take N (number of sampling points of the function f), x_inf and x_sup (the limits of the domain of f)

### OUTPUT

- save a file with N row and two columns representing the value of x and the respective f(x) suing the information from the input

- print on terminal the value (with 16 decimal digits) of the following integral, calculated using trapezoidal rule:

  $I = \int_{0}^{\pi/2} f(x)\,dx$

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>

/* Costante pi portabile */
#define MY_PI 3.14159265358979323846

/* f(x) = e^x * cos(x) */
static inline double f(double x) {
    return exp(x) * cos(x);
}

/* Regola dei trapezi su [a,b] con M sottointervalli (M >= 1) */
double trapezoidal_integral(double (*func)(double), double a, double b, long long M) {
    if (M < 1) M = 1;
    const double h = (b - a) / (double)M;

    double sum = 0.5 * (func(a) + func(b));
    for (long long k = 1; k < M; ++k) {
        double xk = a + h * (double)k;
        sum += func(xk);
    }
    return h * sum;
}

/* --- Utility: trim spazi in-place --- */
static void trim_spaces(char *s) {
    // remove leading spaces
    char *p = s;
    while (isspace((unsigned char)*p)) p++;
    if (p != s) memmove(s, p, strlen(p) + 1);

    // remove trailing spaces
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) {
        s[n - 1] = '\0';
        n--;
    }
}

/* --- Utility: remove ALL spaces (in-place) --- */
static void remove_all_spaces(char *s) {
    char *dst = s;
    for (char *src = s; *src; ++src) {
        if (!isspace((unsigned char)*src)) {
            *dst++ = *src;
        }
    }
    *dst = '\0';
}

/*
  Parse numeri tipo:
  - "1.23"
  - "pi", "-pi"
  - "pi/2", "-pi/2"
  - "2*pi", "3*pi/4", "1.5*pi"
  - "pi*2" (anche questa)
  - "2*pi/3"
  Regole:
  - ignora spazi
  - supporta un solo 'pi' e una sola divisione '/'
*/
int parse_math_expr(const char *input, double *out) {
    char buf[256];
    if (strlen(input) >= sizeof(buf)) return 0;
    strcpy(buf, input);
    trim_spaces(buf);
    remove_all_spaces(buf);

    if (buf[0] == '\0') return 0;

    // Gestisco eventuale divisione "/"
    char *slash = strchr(buf, '/');
    char left[256], right[256];
    if (slash) {
        *slash = '\0';
        strncpy(left, buf, sizeof(left));
        left[sizeof(left)-1] = '\0';
        strncpy(right, slash + 1, sizeof(right));
        right[sizeof(right)-1] = '\0';
    } else {
        strncpy(left, buf, sizeof(left));
        left[sizeof(left)-1] = '\0';
        right[0] = '\0';
    }

    // Funzione interna: valuta "term" che può contenere pi e moltiplicazioni
    // Esempi: "pi", "-pi", "2*pi", "pi*2", "1.5*pi", "-3*pi"
    // Oppure solo numero: "2.7"
    auto double eval_term(const char *term, int *ok) {
        *ok = 1;
        if (term[0] == '\0') { *ok = 0; return 0.0; }

        // Cerca "pi"
        const char *pi_pos = strstr(term, "pi");
        if (!pi_pos) {
            // solo numero
            char *endp = NULL;
            errno = 0;
            double v = strtod(term, &endp);
            if (errno != 0 || endp == term || *endp != '\0') { *ok = 0; return 0.0; }
            return v;
        }

        // Se c'è pi, sostituiamo: valore = coeff * MY_PI
        // term può essere:
        //   "pi" / "-pi"
        //   "<coeff>*pi"  o "<coeff>pi" (accettiamo anche senza '*')
        //   "pi*<coeff>"
        //   "<coeff>*pi*<coeff2>" -> non supportato (teniamolo semplice)
        char tmp[256];
        strncpy(tmp, term, sizeof(tmp));
        tmp[sizeof(tmp)-1] = '\0';

        // Caso "pi" o "-pi"
        if (strcmp(tmp, "pi") == 0) return MY_PI;
        if (strcmp(tmp, "-pi") == 0) return -MY_PI;

        // Caso "pi*X"
        if (strncmp(tmp, "pi*", 3) == 0) {
            const char *rhs = tmp + 3;
            char *endp = NULL;
            errno = 0;
            double c = strtod(rhs, &endp);
            if (errno != 0 || endp == rhs || *endp != '\0') { *ok = 0; return 0.0; }
            return MY_PI * c;
        }
        if (strncmp(tmp, "-pi*", 4) == 0) {
            const char *rhs = tmp + 4;
            char *endp = NULL;
            errno = 0;
            double c = strtod(rhs, &endp);
            if (errno != 0 || endp == rhs || *endp != '\0') { *ok = 0; return 0.0; }
            return -MY_PI * c;
        }

        // Caso "X*pi" o "Xpi" o "X*pi" con segno
        // Rimuoviamo la parte "pi" e optional '*'
        char *p = strstr(tmp, "pi");
        if (!p) { *ok = 0; return 0.0; }
        *p = '\0'; // tronca prima di pi

        // ora tmp contiene qualcosa come "2*" oppure "2" oppure "-3*" ecc.
        size_t len = strlen(tmp);
        if (len > 0 && tmp[len-1] == '*') tmp[len-1] = '\0';

        // se rimane vuoto, coeff = 1
        if (tmp[0] == '\0') return MY_PI;

        char *endp = NULL;
        errno = 0;
        double c = strtod(tmp, &endp);
        if (errno != 0 || endp == tmp || *endp != '\0') { *ok = 0; return 0.0; }
        return c * MY_PI;
    }

    int ok1 = 0, ok2 = 1;
    double v1 = eval_term(left, &ok1);
    if (!ok1) return 0;

    if (slash) {
        double v2 = eval_term(right, &ok2);
        if (!ok2 || v2 == 0.0) return 0;
        *out = v1 / v2;
        return 1;
    } else {
        *out = v1;
        return 1;
    }
}

/* Legge una linea da stdin in modo sicuro */
static int read_line(const char *prompt, char *buf, size_t buflen) {
    printf("%s", prompt);
    fflush(stdout);
    if (!fgets(buf, (int)buflen, stdin)) return 0;
    // rimuovi newline
    size_t n = strlen(buf);
    if (n > 0 && buf[n-1] == '\n') buf[n-1] = '\0';
    return 1;
}

int main(void) {
    char line[256];

    printf("Task04 interactive mode\n");
    printf("Function: f(x) = exp(x) * cos(x)\n\n");

    // N
    long long N = 0;
    while (1) {
        if (!read_line("Enter N (number of sampling points, >= 2): ", line, sizeof(line))) return EXIT_FAILURE;
        char *endp = NULL;
        errno = 0;
        N = strtoll(line, &endp, 10);
        if (errno == 0 && endp != line && *endp == '\0' && N >= 2) break;
        printf("Invalid N. Try again.\n");
    }

    // x_inf
    double x_inf = 0.0;
    while (1) {
        if (!read_line("Enter x_inf (e.g. -1, 0, -pi/2, -3*pi/4): ", line, sizeof(line))) return EXIT_FAILURE;
        if (parse_math_expr(line, &x_inf)) break;
        printf("Invalid x_inf. Try again (allowed: numbers, pi, pi/2, 2*pi, 3*pi/4, ...).\n");
    }

    // x_sup
    double x_sup = 0.0;
    while (1) {
        if (!read_line("Enter x_sup (e.g. 2, pi/2, pi, 3*pi/2): ", line, sizeof(line))) return EXIT_FAILURE;
        if (parse_math_expr(line, &x_sup) && x_sup > x_inf) break;
        printf("Invalid x_sup (must be > x_inf). Try again.\n");
    }

    // output filename
    char outname[256];
    while (1) {
        if (!read_line("Enter output filename (e.g. fx_data.txt): ", outname, sizeof(outname))) return EXIT_FAILURE;
        trim_spaces(outname);
        if (outname[0] != '\0') break;
        printf("Invalid filename. Try again.\n");
    }

    // M for integral
    long long M_for_integral = 1000000LL;
    if (!read_line("Enter M for trapezoidal integral on [0, pi/2] (press ENTER for default 1000000): ", line, sizeof(line)))
        return EXIT_FAILURE;
    trim_spaces(line);
    if (line[0] != '\0') {
        char *endp = NULL;
        errno = 0;
        long long tmp = strtoll(line, &endp, 10);
        if (errno == 0 && endp != line && *endp == '\0' && tmp >= 1) {
            M_for_integral = tmp;
        } else {
            printf("Invalid M, using default 1000000.\n");
        }
    }

    // Write file
    FILE *fout = fopen(outname, "w");
    if (!fout) {
        perror("Error opening output file");
        return EXIT_FAILURE;
    }

    const double h_sample = (x_sup - x_inf) / (double)(N - 1);
    for (long long i = 0; i < N; ++i) {
        double x = x_inf + h_sample * (double)i;
        double fx = f(x);
        fprintf(fout, "%.16f %.16f\n", x, fx);
    }
    fclose(fout);

    // Integral on [0, pi/2] (task requirement)
    const double a = 0.0;
    const double b = 0.5 * MY_PI;
    const double I_num = trapezoidal_integral(f, a, b, M_for_integral);
    const double I_true = (exp(0.5 * MY_PI) - 1.0) / 2.0;
    const double rel_err = I_num / I_true - 1.0;

    printf("\nSaved output file: %s\n", outname);
    printf("Integral on [0, pi/2] using trapezoidal rule:\n");
    printf("Integrale numerico I  = %.16f\n", I_num);
    printf("Integrale esatto  I* = %.16f\n", I_true);
    printf("Errore relativo      = %.16e\n", rel_err);

    return EXIT_SUCCESS;
}

}

```

### Answers to the questions of Task 04

**1) Printed on the terminal:**
```
Integrale numerico I  = 1.9052386904823926

Integrale esatto  I* = 1.9052386904826757

Errore relativo      = -1.4854784069484595e-13
```
**2) I tried to use the long double instead of double, but nothing changed. So I decided to calculate che integral value with a different method, using Simpson (error∼O($h^4$)) method instead the trapezoidal one(error∼O($h^2$)):**
```
- From terminal, with N sublevel = 1000:
  
Integrale numerico I  = 1.905238690482418066575

Integrale esatto  I* = 1.905238690482675827835

Errore relativo      = -1.352908128408880639171e-13
```
- From terminal, with N sublevel = 10000:
```
Integrale numerico I  = 1.905238690482675798670

Integrale esatto  I* = 1.905238690482675827835

Errore relativo      = -1.528725063204561251951e-17
 ```  
- From terminal, with N sublevel = 1000000:
```  
Integrale numerico I  = 1.905238690482675828161

Integrale esatto  I* = 1.905238690482675827835

Errore relativo      = 2.168404344971008868015e-19 minimum error with N=10^6
```
**4) Using interpreted language:**
```   
Integrale dal file (I4)   = 1.9032135801879795

Integrale teorico (Itrue) = 1.9052386904826757

Errore assoluto |I - I4|  = 2.025e-03
```
