#include <time.h>
#include "fq_poly.h"

#include "../nmod_multi_poly.h"
#include "../fmpz_mod_multi_poly.h"
#include "../fq_multi_poly.h"


#include "../ulong_extras/point_fixe.c"
#include "../nmod_poly/fft_bluestein2.c"
#include "../nmod_multi_poly/mod_fermat.c"
#include "../nmod_multi_poly/fft.c"
#include "../nmod_multi_poly/evaluate_nmod_vec_fast.c"
#include "../nmod_multi_poly/multimodular_primes.c"
#include "../nmod_multi_poly/multimodular.c"
#include "../nmod_multi_poly/init.c"
#include "../nmod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/init.c"
#include "../fmpz_mod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/multimodular.c"
#include "../fq_multi_poly/init.c"
#include "../fq_multi_poly/clear.c"
#include "../fq_poly/interpolate_fq_vec_fast.c"
#include "../fq_multi_poly/multimodular.c"
#include "../fq_poly/compose_mod_kedlaya_umans.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o G-fq_poly_compose_mod_kedlaya_umans g-fq_poly_compose_mod_kedlaya_umans.c -lflint -lmpfr -lgmp -lpthread



int main() {

	clock_t c1, c2;
	double c_temps, c_temps2;
	slong n, lim = 10, d, m, degree = 2;
	fq_ctx_t ctx;
	fq_t one;
	fq_poly_t f, g, h, fg, fg2;
	fmpz_t p, big;
	FILE* F;


	flint_rand_t seed;
	flint_randinit(seed);

	fmpz_init(p);
	fmpz_set_ui(p, 65537);

	fmpz_init_set_ui(big, 2);
	fmpz_pow_ui(big, big, 100);
	fmpz_add_ui(big,big,277);

	fmpz_set(p, big);


	fq_ctx_init(ctx, p, degree, "a");

	fq_init(one,ctx);
	fq_one(one,ctx);

	fq_poly_init(f, ctx);
	fq_poly_init(g, ctx);
	fq_poly_init(h, ctx);
	fq_poly_init(fg, ctx);
	fq_poly_init(fg2, ctx);

	F = fopen("../../../Gnuplot/plot_fq_poly_compose_mod_kedlaya_umans", "w");

	flint_fprintf(F,"# Caractéristique : ");
	fmpz_fprint(F, p);
	flint_fprintf(F,"\n# Degré de l'extension : %d\n", degree);
	flint_fprintf(F,"n\t\tK-U\t\t\t\tFLINT\t\t\td\t\tm\n");

	for(n = 5 ; n <= lim ; n++) {
		printf("n : %ld\n", n);

		d = n_sqrt(n)+1;
		m = n_clog(n,d);

		flint_randinit(seed);
		fq_poly_randtest(f, seed, n-1, ctx);
		flint_randinit(seed);
		fq_poly_randtest(g, seed, n-1, ctx);
		flint_randinit(seed);
		fq_poly_randtest(h, seed, n, ctx);
		// Coefficient dominant à 1
		fq_poly_set_coeff(h, n-1, one, ctx);


		c1 = clock();
		fq_poly_compose_mod_kedlaya_umans(fg,f,g,h,d,ctx);
		c2 = clock();

		c_temps = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;


		c1 = clock();
		fq_poly_compose_mod(fg2, f, g, h, ctx);
		c2 = clock();

		c_temps2 = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;

		flint_fprintf(F,"%d\t\t%f\t\t%f\t\t%d\t\t%d\n",n,c_temps,c_temps2,d,m);

		if(fq_poly_equal(fg,fg2,ctx) == 0) {
			printf("Problème de validité !\n");
			exit(1);
		}

	}

	fclose(F);

	flint_randclear(seed);


	fmpz_clear(p);
	fmpz_clear(big);
	fq_clear(one,ctx);
	fq_poly_clear(f,ctx);
	fq_poly_clear(g,ctx);
	fq_poly_clear(h,ctx);
	fq_poly_clear(fg,ctx);
	fq_poly_clear(fg2,ctx);
	fq_ctx_clear(ctx);


	return 0;
}


