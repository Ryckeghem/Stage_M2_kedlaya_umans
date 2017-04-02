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



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-fmpz_mod_multi_poly_multimodular t-fmpz_mod_multi_poly_multimodular.c -lflint -lmpfr -lgmp -lpthread



int main() {

	fmpz_t p;
	slong i, d, m, N, dm, pp;
	fmpz_mod_multi_poly_t f;
	fmpz *vect_alpha, *vect_res;

	fmpz_init_set_ui(p, 2003);
	pp = 2003;
	d = 10;
	m = 2;
	dm = n_pow(d,m);
	N = dm*m*d;

	fmpz_mod_multi_poly_init(f, dm, *p, d, m);
	_fmpz_vec_zero(f->poly, dm);
	fmpz_set_ui(f->poly,0);
	fmpz_set_ui(f->poly +1,0);
	// f = x_2
	fmpz_set_ui(f->poly +10,1);

	vect_alpha = _fmpz_vec_init(N*m);
	for(i = 0 ; i < N*m ; i++) {
		vect_alpha[i] = i%pp;
	}

	vect_res = _fmpz_vec_init(N);
	fmpz_mod_multi_poly_multimodular(vect_res, N, f, vect_alpha);

	_fmpz_vec_print(vect_res,N);
	flint_printf("\n");


	fmpz_clear(p);
	fmpz_mod_multi_poly_clear(f);
	_fmpz_vec_clear(vect_alpha, N*m);
	_fmpz_vec_clear(vect_res, N);

	return 0;
}
