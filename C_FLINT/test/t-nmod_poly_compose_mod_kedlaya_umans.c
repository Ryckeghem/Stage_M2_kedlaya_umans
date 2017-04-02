#include "nmod_poly.h"

#include "../nmod_multi_poly.h"


#include "../ulong_extras/point_fixe.c"
#include "../nmod_poly/fft_bluestein2.c"
#include "../nmod_multi_poly/mod_fermat.c"
#include "../nmod_multi_poly/fft.c"
#include "../nmod_multi_poly/evaluate_nmod_vec_fast.c"
#include "../nmod_multi_poly/multimodular_primes.c"
#include "../nmod_multi_poly/multimodular.c"
#include "../nmod_multi_poly/init.c"
#include "../nmod_multi_poly/clear.c"
#include "../nmod_poly/compose_mod_kedlaya_umans.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-nmod_poly_compose_mod_kedlaya_umans t-nmod_poly_compose_mod_kedlaya_umans.c -lflint -lmpfr -lgmp -lpthread



int main() {

	slong d = 3;
	mp_limb_t p, llbig, lbig, big, big_prime;
	nmod_poly_t f, g, h, fg, fg2;

	llbig = 65535ULL;
	lbig = 4294967295ULL;
	big = 18446744073709551556ULL;
	big_prime = 18446744073709551557ULL;
	p = big_prime;

	nmod_poly_init(f, p);
	nmod_poly_init(g, p);
	nmod_poly_init(h, p);
	nmod_poly_init(fg, p);
	nmod_poly_init(fg2, p);

	nmod_poly_set_coeff_ui(f, 0, 9);
	nmod_poly_set_coeff_ui(f, 1, big);
	nmod_poly_set_coeff_ui(f, 2, 2);
	nmod_poly_set_coeff_ui(f, 3, 3);
	nmod_poly_set_coeff_ui(f, 4, lbig);
	nmod_poly_set_coeff_ui(f, 5, llbig);
	nmod_poly_set_coeff_ui(f, 6, 6);
	nmod_poly_set_coeff_ui(f, 7, 7);

	nmod_poly_set_coeff_ui(g, 0, 9);
	nmod_poly_set_coeff_ui(g, 1, lbig);
	nmod_poly_set_coeff_ui(g, 2, 2);
	nmod_poly_set_coeff_ui(g, 3, 3);
	nmod_poly_set_coeff_ui(g, 4, llbig);
	nmod_poly_set_coeff_ui(g, 5, big);
	nmod_poly_set_coeff_ui(g, 6, 6);
	nmod_poly_set_coeff_ui(g, 7, 7);

	nmod_poly_set_coeff_ui(h, 0, 9);
	nmod_poly_set_coeff_ui(h, 1, p-1);
	nmod_poly_set_coeff_ui(h, 2, 2);
	nmod_poly_set_coeff_ui(h, 3, 3);
	nmod_poly_set_coeff_ui(h, 4, llbig);
	nmod_poly_set_coeff_ui(h, 5, big);
	nmod_poly_set_coeff_ui(h, 6, 6);
	nmod_poly_set_coeff_ui(h, 7, 7);
	nmod_poly_set_coeff_ui(h, 8, 1);

	flint_printf("f : ");	
	nmod_poly_print(f);
	flint_printf("\ng : ");
	nmod_poly_print(g);
	flint_printf("\nh : ");
	nmod_poly_print(h);
	flint_printf("\n");

	nmod_poly_compose_mod_kedlaya_umans(fg,f,g,h,d);

	printf("Res\n");

	nmod_poly_print(fg);
	printf("\n");

	printf("Res exact\n");
	nmod_poly_compose_mod(fg2, f, g, h);

	nmod_poly_print(fg2);
	flint_printf("\n");

	flint_printf("Booléen d'égalité : %d\n",nmod_poly_equal(fg,fg2));


	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(h);
	nmod_poly_clear(fg);
	nmod_poly_clear(fg2);

	return 0;
}
