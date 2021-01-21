#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi(double a, double b, int n, int i, double x)
{
}

/* Pierwsza pochodna fi */
double
dfi(double a, double b, int n, int i, double x)
{
}

/* Druga pochodna fi */
double
d2fi(double a, double b, int n, int i, double x)
{
}

/* Trzecia pochodna fi */
double
d3fi(double a, double b, int n, int i, double x)
{
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
}

void
make_spl(points_t * pts, spline_t * spl)
{
}
