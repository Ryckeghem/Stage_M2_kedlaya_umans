#include <time.h>
#include "fq_nmod_poly.h"

#include "../nmod_multi_poly.h"


#include "../nmod_poly/fft_bluestein.c"
#include "../nmod_poly/fft_bluestein2.c"
#include "../nmod_poly/fft_geometrique.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o G-nmod_poly_fft_bluestein2 g-nmod_poly_fft_bluestein2.c -lflint -lmpfr -lgmp -lpthread



int main() {

	clock_t c1, c2;
	double c_temps, c_temps2, c_temps3, c_temps4;
	slong i, lim = 20000, degree;
	mp_limb_t w, p;
	mp_ptr vect, res, res2, res3, res4;
	nmod_poly_t f, g, h;
	flint_rand_t seed;
	n_primes_t iter;
	FILE* F;


	flint_randinit(seed);

	vect = _nmod_vec_init(lim);
	res = _nmod_vec_init(lim);
	res2 = _nmod_vec_init(lim);
	res3 = _nmod_vec_init(lim);
	res4 = _nmod_vec_init(lim);

	for(i = 0 ; i < lim ; i++) {
		vect[i] = i;
	}

	F = fopen("../../../Gnuplot/plot_nmod_poly_fft_bluestein2", "w");

	fprintf(F,"p\t\tFFT Bluestein\tÉvaluation multipoint\tÉvaluation géométrique\tDegré\n");

	n_primes_init(iter);


	for(p = n_primes_next(iter) ; p <= lim ; p = n_primes_next(iter)) {

		nmod_poly_init(f, p);
		nmod_poly_init(g, p);
		nmod_poly_init(h, p);


		flint_randinit(seed);
		nmod_poly_randtest(f, seed, p);
		flint_randinit(seed);

		degree = nmod_poly_degree(f);

		nmod_poly_set(g, f);
		nmod_poly_set(h, f);


		c1 = clock();
		res2[0] = f->coeffs[0];
		nmod_poly_evaluate_nmod_vec(res2+1,f,vect+1,p-1);
		c2 = clock();

		c_temps2 = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;



		w = n_primitive_root_prime(p);

		// Bluestein2
		c1 = clock();

		if(nmod_poly_length(f) == p) {
			i = f->coeffs[0];
			nmod_poly_set_coeff_ui(f, 0, n_addmod(f->coeffs[0], f->coeffs[p-1], p));
			nmod_poly_set_coeff_ui(f, p-1, 0);

			fft_bluestein2(res,f,w);
			res[0] = i;
		}
		else {
			fft_bluestein2(res,f,w);
		}
		c2 = clock();

		c_temps = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;


		// Bluestein
		c1 = clock();

		if(nmod_poly_length(g) == p) {
			i = g->coeffs[0];
			nmod_poly_set_coeff_ui(g, 0, n_addmod(g->coeffs[0], g->coeffs[p-1], p));
			nmod_poly_set_coeff_ui(g, p-1, 0);

			fft_bluestein(res3,g,w);
			res3[0] = i;
		}
		else {
			fft_bluestein(res3,g,w);
		}
		c2 = clock();

		c_temps3 = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;


		// Version géométrique
		c1 = clock();

		if(nmod_poly_length(h) == p) {
			i = h->coeffs[0];
			nmod_poly_set_coeff_ui(h, 0, n_addmod(h->coeffs[0], h->coeffs[p-1], p));
			nmod_poly_set_coeff_ui(h, p-1, 0);

			fft_geometrique(res4,h,w);
			res4[0] = i;
		}
		else {
			fft_geometrique(res4,h,w);
		}
		c2 = clock();

		c_temps4 = ((double)c2 - (double)c1)/CLOCKS_PER_SEC;



		flint_fprintf(F,"%ld\t\t%f\t\t%f\t\t%f\t\t%f\t\t\t\t%ld\n",p,c_temps,c_temps3,c_temps2,c_temps4,degree);

		if((_nmod_vec_equal(res2,res,p) == 0) || (_nmod_vec_equal(res2,res3,p) == 0) || (_nmod_vec_equal(res2,res4,p) == 0)) {
			printf("Problème de validité !\n");
			exit(1);
		}

	}

	fclose(F);

	flint_randclear(seed);
	n_primes_clear(iter);


	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(h);
	_nmod_vec_clear(vect);
	_nmod_vec_clear(res);
	_nmod_vec_clear(res2);
	_nmod_vec_clear(res3);
	_nmod_vec_clear(res4);


	return 0;
}
