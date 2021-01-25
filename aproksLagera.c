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

double liczPoch(double x, int i, int st){

double sum=0;
int it;

if(x==0 )
	sum=0;
else{
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
		
		double tmp = ( (jedynka/silniaK) * iks * silniaI / (silniaK * silniaNK ) ) ;
		
		int gj;
		int c=k;
//		printf("k: %d\n\n", k);
		for( gj=0; gj<st; gj++){
			tmp = tmp*c/x;
			c--;
		}
		sum += tmp;
	}
}
return sum;
}


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
return liczPoch(x, i, 1);
}

/* Druga pochodna fi */
double
d2fi(double a, double b, int n, int i, double x)
{
return liczPoch(x, i, 2);
}

/* Trzecia pochodna fi */
double
d3fi(double a, double b, int n, int i, double x)
{
return liczPoch(x, i, 3);
}

double
xfi(double a, double b, int n, int i, FILE *out)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %d", hi[j] );
	fprintf( out, "] hx=[" );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %g", hx[j] );
	fprintf( out, "]\n" );
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);


	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
	}


	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}


}
