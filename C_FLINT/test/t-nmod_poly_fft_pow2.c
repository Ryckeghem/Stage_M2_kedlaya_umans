#include "nmod_poly.h"


#include "../nmod_poly/fft_pow2.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-nmod_poly_fft_pow2 t-nmod_poly_fft_pow2.c -lflint -lmpfr -lgmp -lpthread



int main() {

	slong i;
	mp_limb_t p, w;
	nmod_poly_t f;
	mp_ptr res, res2, vect;

	// p premier
	// Le cardinal du groupe multiplicatif doit être une puissance de 2
	// Ne marche que pour p = 0, 1, 2, 3, 5, 17, 257, 65537.
	p = 17;

	nmod_poly_init(f, p);
	nmod_poly_set_coeff_ui(f, 0, 0);
	nmod_poly_set_coeff_ui(f, 1, 1);
	nmod_poly_set_coeff_ui(f, 2, 2);
	nmod_poly_set_coeff_ui(f, 3, 1);
	nmod_poly_set_coeff_ui(f, 4, 1);

	w = n_primitive_root_prime(p);

	res = _nmod_vec_init(p);
	nmod_poly_fft_pow2(res, f, w);

	flint_printf("w : %d\n", w);
	_fmpz_vec_print(res, p);
	flint_printf("\n");


	vect = _nmod_vec_init(p);
	for(i = 0 ; i < p ; i++) {
		vect[i] = i;
	}
	res2 = _nmod_vec_init(p);
	nmod_poly_evaluate_nmod_vec(res2, f, vect, p);

	_fmpz_vec_print(res2, p);
	flint_printf("\nBooléen d'égalité : %d\n", _nmod_vec_equal(res,res2,p));

	nmod_poly_clear(f);
	_nmod_vec_clear(res);
	_nmod_vec_clear(res2);
	_nmod_vec_clear(vect);

	return 0;
}
