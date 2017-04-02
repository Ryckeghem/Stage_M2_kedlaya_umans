#include "fmpz_mod_poly.h"

#include "../nmod_multi_poly.h"
#include "../fmpz_mod_multi_poly.h"


#include "../ulong_extras/point_fixe.c"
#include "../nmod_poly/fft_bluestein2.c"
#include "../nmod_multi_poly/clear.c"
#include "../nmod_multi_poly/mod_fermat.c"
#include "../nmod_multi_poly/fft.c"
#include "../nmod_multi_poly/evaluate_nmod_vec_fast.c"
#include "../nmod_multi_poly/multimodular_primes.c"
#include "../fmpz_mod_multi_poly/init.c"
#include "../fmpz_mod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/multimodular.c"
#include "../fmpz_mod_poly/interpolate_fmpz_vec_fast.c"
#include "../fmpz_mod_poly/compose_mod_kedlaya_umans.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-fmpz_mod_poly_compose_mod_kedlaya_umans t-fmpz_mod_poly_compose_mod_kedlaya_umans.c -lflint -lmpfr -lgmp -lpthread



int main() {

	slong d = 3;
	fmpz_t p, big, nb;
	fmpz_mod_poly_t f, g, h, fg, fg2;


	// 2^200 + 235 est premier
	fmpz_init_set_ui(big, 2);
	fmpz_pow_ui(big, big, 200);
	fmpz_add_ui(big,big,235);

	fmpz_init_set_ui(nb,2);
	fmpz_pow_ui(nb, nb, 101);

	fmpz_set(p, big);
	fmpz_mod_poly_init(f, p);
	fmpz_mod_poly_init(g, p);
	fmpz_mod_poly_init(h, p);
	fmpz_mod_poly_init(fg, p);
	fmpz_mod_poly_init(fg2, p);

	fmpz_mod_poly_set_coeff_fmpz(f, 0, nb);
	fmpz_mod_poly_set_coeff_ui(f, 1, 1);
	fmpz_mod_poly_set_coeff_ui(f, 2, 2);
	fmpz_mod_poly_set_coeff_ui(f, 3, 3);
	fmpz_mod_poly_set_coeff_ui(f, 4, 4);
	fmpz_mod_poly_set_coeff_ui(f, 5, 5);
	fmpz_mod_poly_set_coeff_fmpz(f, 6, nb);
	fmpz_mod_poly_set_coeff_ui(f, 7, 7);


	fmpz_mod_poly_set_coeff_ui(g, 0, 9);
	fmpz_mod_poly_set_coeff_ui(g, 1, 1);
	fmpz_mod_poly_set_coeff_ui(g, 2, 2);
	fmpz_mod_poly_set_coeff_fmpz(g, 3, nb);
	fmpz_mod_poly_set_coeff_ui(g, 4, 4);
	fmpz_mod_poly_set_coeff_ui(g, 5, 5);
	fmpz_mod_poly_set_coeff_ui(g, 6, 6);
	fmpz_mod_poly_set_coeff_fmpz(g, 7, nb);


	fmpz_mod_poly_set_coeff_ui(h, 0, 9);
	fmpz_mod_poly_set_coeff_si(h, 1, -1);
	fmpz_mod_poly_set_coeff_ui(h, 2, 2);
	fmpz_mod_poly_set_coeff_fmpz(h, 3, nb);
	fmpz_mod_poly_set_coeff_fmpz(h, 4, nb);
	fmpz_mod_poly_set_coeff_ui(h, 5, 5);
	fmpz_mod_poly_set_coeff_ui(h, 6, 6);
	fmpz_mod_poly_set_coeff_ui(h, 7, 7);
	fmpz_mod_poly_set_coeff_ui(h, 8, 1);


	flint_printf("f : ");	
	fmpz_mod_poly_print(f);
	flint_printf("\ng : ");
	fmpz_mod_poly_print(g);
	flint_printf("\nh : ");
	fmpz_mod_poly_print(h);
	flint_printf("\n");

	fmpz_mod_poly_compose_mod_kedlaya_umans(fg,f,g,h,d);

	printf("Res\n");

	fmpz_mod_poly_print(fg);
	printf("\n");

	printf("Res exact\n");
	fmpz_mod_poly_compose_mod(fg2, f, g, h);

	fmpz_mod_poly_print(fg2);
	flint_printf("\n");

	flint_printf("Booléen d'égalité : %d\n",fmpz_mod_poly_equal(fg,fg2));


	fmpz_clear(p);
	fmpz_clear(big);
	fmpz_clear(nb);
	fmpz_mod_poly_clear(f);
	fmpz_mod_poly_clear(g);
	fmpz_mod_poly_clear(h);
	fmpz_mod_poly_clear(fg);
	fmpz_mod_poly_clear(fg2);


	return 0;
}
