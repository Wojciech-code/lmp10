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

//Funkcja obliczająca wartość x w funkcji podstawowej danego stopnia z bazy
double
fi(double a, double b, int n, int i, double x)
{

double sum=1;
int it;

double silniaI=1;
for(it=1; it<=i; it++){
	silniaI=silniaI*it;
}

int k;
for(k=1; k<=i; k++){
	double silniaK=1;
	double jedynka=1;
	double iks=1;
	int it;
	for(it=1; it<=k; it++){
		silniaK=silniaK*it;
		jedynka=jedynka*(-1);
		iks=iks*x;
	}

	double silniaNK= 1;
	int j;
	for( j=1; j<=(i - k); j++){
		silniaNK = silniaNK * j; 
	}

	sum+=(jedynka/silniaK) * iks * silniaI / (silniaK * silniaNK );
}

return sum;
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