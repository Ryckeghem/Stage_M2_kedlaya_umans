#include <time.h>
#include "fmpz_mod_poly.h"

#include "../nmod_multi_poly.h"
#include "../fmpz_mod_multi_poly.h"


#include "../ulong_extras/point_fixe.c"
#include "../nmod_poly/fft_bluestein2.c"
#include "../nmod_multi_poly/mod_fermat.c"
#include "../nmod_multi_poly/fft.c"
#include "../nmod_multi_poly/evaluate_nmod_vec_fast.c"
#include "../nmod_multi_poly/multimodular_primes.c"
#include "../fmpz_mod_multi_poly/multimodular.c"
#include "../nmod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/init.c"
#include "../fmpz_mod_multi_poly/clear.c"
#include "../fmpz_mod_poly/interpolate_fmpz_vec_fast.c"
#include "../fmpz_mod_poly/compose_mod_kedlaya_umans.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o G-fmpz_mod_poly_compose_mod_kedlaya_umans g-fmpz_mod_poly_compose_mod_kedlaya_umans.c -lflint -lmpfr -lgmp -lpthread



int main() {

	clock_t c1, c2;
	double c_temps, c_temps2;
	slong n, lim = 20, d, m, m_max = 2;
	fmpz_mod_poly_t f, g, h, fg, fg2;
	fmpz_t p, big, x, y;
	flint_rand_t seed;
	FILE* F;


	flint_randinit(seed);

	fmpz_init_set_ui(big, 2);
	fmpz_pow_ui(big, big, 200);
	fmpz_add_ui(big,big,235);

	fmpz_init(p);
	fmpz_set(p, big);


	fmpz_mod_poly_init(f, p);
	fmpz_mod_poly_init(g, p);
	fmpz_mod_poly_init(h, p);
	fmpz_mod_poly_init(fg, p);
	fmpz_mod_poly_init(fg2, p);

	F = fopen("../../../Gnuplot/plot_fmpz_mod_poly_compose_mod_kedlaya_umans", "w");

	fprintf(F,"# Caractéristique : ");
	fmpz_fprint(F, p);
	fprintf(F,"\nn\t\t\t\t\tK-U\t\t\t\tFLINT\t\t\td\t\tm\n");

	fmpz_init(x);
	fmpz_init(y);

	for(n = 5 ; n <= lim ; n++) {
		printf("n : %ld\n", n);

		if(m_max == 2) {
			d = n_sqrt(n)+1;
			m = n_clog(n,d);
		}
		else {
			fmpz_set_ui(y,n);
			fmpz_root(x, y, m_max);
			d = fmpz_get_ui(x)+1;
			if(n_pow(d-1,m_max) == n) d--;
			m = n_clog(n,d);
		}


		flint_randinit(seed);
		fmpz_mod_poly_randtest(f, seed, n-1);
		flint_randinit(seed);
		fmpz_mod_poly_randtest(g, seed, n-1);
		flint_randinit(seed);
		fmpz_mod_poly_randtest(h, seed, n);
		// Coefficient dominant à 1
		fmpz_mod_poly_set_coeff_ui(h, n-1, 1);


		c1 = clock();
		fmpz_mod_poly_compose_mod_kedlaya_umans(fg,f,g,h,d);
		c2 = clock();

		c_temps = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;


		c1 = clock();
		fmpz_mod_poly_compose_mod(fg2, f, g, h);
		c2 = clock();

		c_temps2 = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;

		flint_fprintf(F,"%d\t\t\t\t\t%f\t\t%f\t\t%d\t\t%d\n",n,c_temps,c_temps2,d,m);

		fflush(F);

		if(fmpz_mod_poly_equal(fg,fg2) == 0) {
			printf("Problème de validité !\n");
			exit(1);
		}

	}

	fclose(F);

	flint_randclear(seed);


	fmpz_clear(x);
	fmpz_clear(y);
	fmpz_clear(p);
	fmpz_clear(big);
	fmpz_mod_poly_clear(f);
	fmpz_mod_poly_clear(g);
	fmpz_mod_poly_clear(h);
	fmpz_mod_poly_clear(fg);
	fmpz_mod_poly_clear(fg2);


	return 0;
}


