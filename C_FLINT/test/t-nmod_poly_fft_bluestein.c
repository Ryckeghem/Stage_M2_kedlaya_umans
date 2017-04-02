#include "fq_nmod_poly.h"


#include "../nmod_poly/fft_bluestein.c"
#include "../nmod_poly/fft_bluestein2.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-nmod_poly_fft_bluestein t-nmod_poly_fft_bluestein.c -lflint -lmpfr -lgmp -lpthread



int main() {

	flint_rand_t seed;

	slong i;
	mp_limb_t p = 7, w = n_primitive_root_prime(p);
	nmod_poly_t g;
	mp_ptr res, res2, res3, val;


	nmod_poly_init(g, p);
	flint_randinit(seed);
	// degree <= p-2
	// length <= p-1
	nmod_poly_randtest(g, seed, p-1);


	flint_printf("w : %d\n",w);
	printf("g : ");
	nmod_poly_print(g);
	printf("\n");


	res = _nmod_vec_init(p);
	res3 = _nmod_vec_init(p);
	// Choix de la fonction
	fft_bluestein(res,g,w);
	fft_bluestein2(res3,g,w);


	res2 = _nmod_vec_init(p);
	val = _nmod_vec_init(p);
	for(i = 0 ; i < p ; i++) {
		val[i] = i;
	}
	nmod_poly_evaluate_nmod_vec(res2, g, val, p);


	printf("Naif :       ");
	_fmpz_vec_print(res2,p);
	printf("\n");
	printf("Bluestein :  ");
	_fmpz_vec_print(res,p);
	printf("\n");
	flint_printf("Booléen d'égalité : %d\n", _nmod_vec_equal(res,res2,p));
	printf("Bluestein2 : ");
	_fmpz_vec_print(res3,p);
	printf("\n");
	flint_printf("Booléen d'égalité : %d\n", _nmod_vec_equal(res3,res2,p));

	nmod_poly_clear(g);
	_nmod_vec_clear(res);
	_nmod_vec_clear(res2);
	_nmod_vec_clear(res3);
	_nmod_vec_clear(val);

	return 0;
}

