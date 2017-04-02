#include "fmpz_mod_poly.h"


#include "../fmpz_mod_poly/interpolate_fmpz_vec_fast.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-fmpz_mod_poly_interpolate_fmpz_vec_fast t-fmpz_mod_poly_interpolate_fmpz_vec_fast.c -lflint -lmpfr -lgmp -lpthread



int main() {

	slong n = 7;
	fmpz_t mod;
	fmpz* xs = _fmpz_vec_init(n);
	fmpz* ys = _fmpz_vec_init(n);
	fmpz_mod_poly_t res;
	fmpz_poly_t res2;


	fmpz_init(mod);
	fmpz_set_ui(mod, 11);


	fmpz_set_ui(xs, 0);
	fmpz_set_ui(xs+1, 1);
	fmpz_set_ui(xs+2, 2);
	fmpz_set_ui(xs+3, 3);
	fmpz_set_ui(xs+4, 4);
	fmpz_set_ui(xs+5, 5);
	fmpz_set_ui(xs+6, 6);

	fmpz_set_si(ys, 5);
	fmpz_set_si(ys+1, 1);
	fmpz_set_si(ys+2, 4);
	fmpz_set_si(ys+3, 8);
	fmpz_set_si(ys+4, 4);
	fmpz_set_si(ys+5, 1);
	fmpz_set_si(ys+6, 5);

	printf("xs :\n");
	_fmpz_vec_print(xs,n);
	printf("\n");	

	printf("ys :\n");
	_fmpz_vec_print(ys,n);
	printf("\n");

	fmpz_mod_poly_init(res,mod);
	fmpz_mod_poly_interpolate_fmpz_vec_fast(res, xs, ys, n, mod);

	fmpz_poly_init(res2);
	fmpz_poly_interpolate_fmpz_vec(res2, xs, ys, n);

	printf("Nouvelle interpolation (f(xs) = ys) :\n");
	fmpz_mod_poly_print(res);
	printf("\n");

	fmpz* zs = _fmpz_vec_init(n);
	fmpz_poly_evaluate_fmpz_vec(zs, res2, xs, n);

	fmpz* as = _fmpz_vec_init(n);
	fmpz_mod_poly_evaluate_fmpz_vec(as, res, xs, n);

	printf("f(xs) :\n");
	_fmpz_vec_print(as,n);
	printf("\n");

	printf("FLINT :\n");
	fmpz_poly_print(res2);
	printf("\n");

	printf("f(xs) :\n");
	_fmpz_vec_print(zs,n);
	printf("\n");


	fmpz_clear(mod);
	_fmpz_vec_clear(xs,n);
	_fmpz_vec_clear(ys,n);

	return 0;
}
