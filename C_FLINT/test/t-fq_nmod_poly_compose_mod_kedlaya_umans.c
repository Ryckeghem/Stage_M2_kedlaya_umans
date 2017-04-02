#include "fq_poly.h"

#include "../nmod_multi_poly.h"
#include "../fmpz_mod_multi_poly.h"
#include "../fq_nmod_multi_poly.h"
#include "../fq_multi_poly.h"


#include "../ulong_extras/point_fixe.c"
#include "../nmod_poly/fft_bluestein2.c"
#include "../nmod_multi_poly/mod_fermat.c"
#include "../nmod_multi_poly/fft.c"
#include "../nmod_multi_poly/evaluate_nmod_vec_fast.c"
#include "../nmod_multi_poly/multimodular_primes.c"
#include "../nmod_multi_poly/init.c"
#include "../nmod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/init.c"
#include "../fmpz_mod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/multimodular.c"
#include "../fq_nmod_multi_poly/init.c"
#include "../fq_nmod_multi_poly/clear.c"
#include "../fq_multi_poly/init.c"
#include "../fq_multi_poly/clear.c"
#include "../fq_nmod_poly/interpolate_fq_nmod_vec_fast.c"
#include "../fq_multi_poly/multimodular.c"
#include "../fq_nmod_multi_poly/multimodular.c"
#include "../fq_nmod_poly/compose_mod_kedlaya_umans.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-fq_nmod_poly_compose_mod_kedlaya_umans t-fq_nmod_poly_compose_mod_kedlaya_umans.c -lflint -lmpfr -lgmp -lpthread



// Le test montre qu'en fonction du degré de l'extension, le résultat est bon ou faux ...
// Degré 11 bon, degré 10 faux
int main() {

	slong d = 3;
	mp_limb_t llbig, lbig, big, big_prime, pp;
	fq_nmod_poly_t f, g, h, fg, fg2;
	fq_nmod_ctx_t ctx;
	fmpz_t p;
	fq_nmod_t f0, f1, f2, f3, f4, f5, f6, f7, fp, fa;


	fmpz_init(p);
	fmpz_set_ui(p, 2);
	fq_nmod_ctx_init(ctx, p, 11, "a");
//	fq_nmod_poly_init(poly, ctx);
	pp = fmpz_get_ui(p);



	llbig = 65535ULL;
	lbig = 4294967295ULL;
//	big = 18446744073709551615ULL;
	big = 18446744073709551556ULL;
	big_prime = 18446744073709551557ULL;
//	p = big_prime;
//	p = 307/*997*/;


	fq_nmod_init(f0,ctx);
	fq_nmod_init(f1,ctx);
	fq_nmod_init(f2,ctx);
	fq_nmod_init(f3,ctx);
	fq_nmod_init(f4,ctx);
	fq_nmod_init(f5,ctx);
	fq_nmod_init(f6,ctx);
	fq_nmod_init(f7,ctx);
	fq_nmod_init(fp,ctx);
	fq_nmod_init(fa,ctx);

	fq_nmod_set_ui(f0,9,ctx);
	fq_nmod_set_ui(f1,1,ctx);
	fq_nmod_set_ui(f2,2,ctx);
	fq_nmod_set_ui(f3,3,ctx);
	fq_nmod_set_ui(f4,4,ctx);
	fq_nmod_set_ui(f5,5,ctx);
	fq_nmod_set_ui(f6,6,ctx);
	fq_nmod_set_ui(f7,7,ctx);
	fq_nmod_set_ui(fp,pp-1,ctx);
	fq_nmod_gen(fa,ctx);


	fq_nmod_poly_init(f, ctx);
	fq_nmod_poly_init(g, ctx);
	fq_nmod_poly_init(h, ctx);
	fq_nmod_poly_init(fg, ctx);
	fq_nmod_poly_init(fg2, ctx);

	fq_nmod_poly_set_coeff(f, 0, f0, ctx);
	fq_nmod_poly_set_coeff(f, 1, fa, ctx);
	fq_nmod_poly_set_coeff(f, 2, f2, ctx);
	fq_nmod_poly_set_coeff(f, 3, f3, ctx);
	fq_nmod_poly_set_coeff(f, 4, f4, ctx);
	fq_nmod_poly_set_coeff(f, 5, f5, ctx);
	fq_nmod_poly_set_coeff(f, 6, f6, ctx);
	fq_nmod_poly_set_coeff(f, 7, f7, ctx);

/*
	fq_nmod_poly_set_coeff_ui(g, 0, 9);
	fq_nmod_poly_set_coeff_ui(g, 1, lbig);
	fq_nmod_poly_set_coeff_ui(g, 2, 2);
	fq_nmod_poly_set_coeff_ui(g, 3, 3);
	fq_nmod_poly_set_coeff_ui(g, 4, llbig);
	fq_nmod_poly_set_coeff_ui(g, 5, big);
	fq_nmod_poly_set_coeff_ui(g, 6, 6);
	fq_nmod_poly_set_coeff_ui(g, 7, 7);

	fq_nmod_poly_set_coeff_ui(h, 0, 9);
	fq_nmod_poly_set_coeff_ui(h, 1, big);
	fq_nmod_poly_set_coeff_ui(h, 2, 2);
	fq_nmod_poly_set_coeff_ui(h, 3, 3);
	fq_nmod_poly_set_coeff_ui(h, 4, llbig);
	fq_nmod_poly_set_coeff_ui(h, 5, big);
	fq_nmod_poly_set_coeff_ui(h, 6, 6);
	fq_nmod_poly_set_coeff_ui(h, 7, 7);
*/
	flint_printf("f : ");	
	fq_nmod_poly_print(f,ctx);
	flint_printf("\ng : ");
	fq_nmod_poly_set_coeff(g, 2, fa, ctx);
	fq_nmod_poly_print(g,ctx);
	flint_printf("\nh : ");
	fq_nmod_poly_set_coeff(h, 1, fa, ctx);
	fq_nmod_poly_set_coeff(h, 8, f1, ctx);
	fq_nmod_poly_print(h,ctx);
	flint_printf("\n");

	fq_nmod_poly_compose_mod_kedlaya_umans(fg,f,g,h,d,ctx);

	printf("Res\n");

	fq_nmod_poly_print(fg,ctx);
	printf("\n");

	printf("Vrai res\n");
	fq_nmod_poly_compose_mod(fg2, f, g, h, ctx);

	fq_nmod_poly_print(fg2,ctx);
	flint_printf("\n");

	flint_printf("Booléen d'égalité : %d\n",fq_nmod_poly_equal(fg,fg2,ctx));


	fmpz_clear(p);
	fq_nmod_poly_clear(f,ctx);
	fq_nmod_poly_clear(g,ctx);
	fq_nmod_poly_clear(h,ctx);
	fq_nmod_poly_clear(fg,ctx);
	fq_nmod_poly_clear(fg2,ctx);
	fq_nmod_clear(f0,ctx);
	fq_nmod_clear(f1,ctx);
	fq_nmod_clear(f2,ctx);
	fq_nmod_clear(f3,ctx);
	fq_nmod_clear(f4,ctx);
	fq_nmod_clear(f5,ctx);
	fq_nmod_clear(f6,ctx);
	fq_nmod_clear(f7,ctx);
	fq_nmod_clear(fp,ctx);
	fq_nmod_clear(fa,ctx);
	fq_nmod_ctx_clear(ctx);

	return 0;
}
