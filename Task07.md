### Task 07

## Interpreted language python:

The function daxpy to test:

```py
# daxpy.py

def daxpy(a, x, y):
    if len(x) != len(y):
        raise ValueError("Vectors must have same length")
    
    d = [0.0] * len(x)
    for i in range(len(x)):
        d[i] = a * x[i] + y[i]
    return d
```

Code for unit test:

```py
import unittest
from daxpy import daxpy

class TestDaxpy(unittest.TestCase):

    def test_small_vector(self):
        a = 3
        x = [1, 2, 3]
        y = [4, 5, 6]
        expected = [a*x[i] + y[i] for i in range(len(x))]
        self.assertEqual(daxpy(a, x, y), expected)

    def test_zero_length(self):
        self.assertEqual(daxpy(3, [], []), [])

    def test_negative_values(self):
        a = -2
        x = [1, -1]
        y = [3, 3]
        expected = [a*x[i] + y[i] for i in range(len(x))]
        self.assertEqual(daxpy(a, x, y), expected)

    def test_size_mismatch(self):
        with self.assertRaises(ValueError):
            daxpy(3, [1,2], [1])

    def test_large(self):
        a = 3
        N = 100000
        x = [1.0]*N
        y = [2.0]*N
        d = daxpy(a, x, y)
        self.assertEqual(d[0], 5.0)
        self.assertEqual(d[N//2], 5.0)
        self.assertEqual(d[-1], 5.0)

if __name__ == '__main__':
    unittest.main()

```
# Output

```
.....
----------------------------------------------------------------------
Ran 5 tests in 0.022s

OK

```
## Compiled language C:

# define first the header file daxpy.h

```c

#ifndef DAXPY_H
#define DAXPY_H

void daxpy(double a, const double *x, const double *y, double *d, int n);

#endif

```

# define now daxpy.c

```c

// daxpy.c
#include "daxpy.h"

void daxpy(double a, const double *x, const double *y, double *d, int n) {
    for (int i = 0; i < n; i++) {
        d[i] = a * x[i] + y[i];
    }
}

```

# test_daxpy.c

```c

// test_daxpy.c
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "daxpy.h"

void test_small() {
    double x[3] = {1, 2, 3};
    double y[3] = {4, 5, 6};
    double d[3];

    daxpy(3.0, x, y, d, 3);

    assert(d[0] == 7.0);
    assert(d[1] == 11.0);
    assert(d[2] == 15.0);
}

void test_zero_length() {
    // n = 0 → la funzione non deve crashare
    daxpy(3.0, NULL, NULL, NULL, 0);
}

void test_negative_values() {
    double x[2] = {1.0, -1.0};
    double y[2] = {3.0, 3.0};
    double d[2];

    daxpy(-2.0, x, y, d, 2);

    assert(d[0] == 1.0);  // -2*1 + 3 = 1
    assert(d[1] == 5.0);  // -2*(-1) + 3 = 5
}

void test_large() {
    int N = 100000;
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));

    assert(x != NULL && y != NULL && d != NULL);

    for (int i = 0; i < N; i++) {
        x[i] = 1.0;
        y[i] = 2.0;
    }

    daxpy(3.0, x, y, d, N);

    assert(d[0] == 5.0);
    assert(d[N/2] == 5.0);
    assert(d[N-1] == 5.0);

    free(x);
    free(y);
    free(d);
}

int main(void) {
    test_small();
    test_zero_length();
    test_negative_values();
    test_large();

    printf("All C daxpy tests passed!\n");
    return 0;
}

```
# We write the makefile

```c
CC = gcc
CFLAGS = -Wall -O2
TARGET = test_daxpy

all: $(TARGET)

$(TARGET): daxpy.o test_daxpy.o
	$(CC) $(CFLAGS) daxpy.o test_daxpy.o -o $(TARGET)

daxpy.o: daxpy.c daxpy.h
	$(CC) $(CFLAGS) -c daxpy.c

test_daxpy.o: test_daxpy.c daxpy.h
	$(CC) $(CFLAGS) -c test_daxpy.c

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o $(TARGET)

```

# Output

```
[root@be0669f7bc73 task07]# make
gcc -Wall -O2 -c daxpy.c
gcc -Wall -O2 -c test_daxpy.c
gcc -Wall -O2 daxpy.o test_daxpy.o -o test_daxpy
[root@be0669f7bc73 task07]# ./test_daxpy 
All C daxpy tests passed!

```
