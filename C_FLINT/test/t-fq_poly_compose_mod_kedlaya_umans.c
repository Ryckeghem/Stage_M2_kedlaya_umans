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
#include "../nmod_poly/compose_mod_kedlaya_umans.c"
#include "../fmpz_mod_multi_poly/init.c"
#include "../fmpz_mod_multi_poly/clear.c"
#include "../fmpz_mod_multi_poly/multimodular.c"
#include "../fq_multi_poly/init.c"
#include "../fq_multi_poly/clear.c"
#include "../fq_poly/interpolate_fq_vec_fast.c"
#include "../fq_multi_poly/multimodular.c"
#include "../fq_poly/compose_mod_kedlaya_umans.c"



// Compilation :
// gcc -I /home/jocelyn/local/include/ -I /home/jocelyn/local/include/flint -L /home/jocelyn/local/lib/ -o T-fq_poly_compose_mod_kedlaya_umans t-fq_poly_compose_mod_kedlaya_umans.c -lflint -lmpfr -lgmp -lpthread



int main() {

	slong d = 3, degree = 4;
	fq_poly_t f, g, h, fg, fg2;
	fq_ctx_t ctx;
	fmpz_t p, p1, big;
	fq_t f0, f1, f2, f3, f4, f5, f6, f7, fp, fa;


	fmpz_init(p);
	fmpz_set_ui(p, 59);

	fmpz_sub_ui(p1,p,1);

	// 2^200 + 235 est premier
	fmpz_init_set_ui(big, 2);
	fmpz_pow_ui(big, big, 200);
	fmpz_add_ui(big,big,235);
	fmpz_set(p, big);

	fq_ctx_init(ctx, p, degree, "a");


	fq_init(f0,ctx);
	fq_init(f1,ctx);
	fq_init(f2,ctx);
	fq_init(f3,ctx);
	fq_init(f4,ctx);
	fq_init(f5,ctx);
	fq_init(f6,ctx);
	fq_init(f7,ctx);
	fq_init(fp,ctx);
	fq_init(fa,ctx);

	fq_set_ui(f0,9,ctx);
	fq_set_ui(f1,1,ctx);
	fq_set_ui(f2,2,ctx);
	fq_set_ui(f3,3,ctx);
	fq_set_ui(f4,4,ctx);
	fq_set_ui(f5,5,ctx);
	fq_set_ui(f6,6,ctx);
	fq_set_ui(f7,7,ctx);
	fq_set_fmpz(fp,p1,ctx);
	fq_gen(fa,ctx);


	fq_poly_init(f, ctx);
	fq_poly_init(g, ctx);
	fq_poly_init(h, ctx);
	fq_poly_init(fg, ctx);
	fq_poly_init(fg2, ctx);

	fq_poly_set_coeff(f, 0, f0, ctx);
	fq_poly_set_coeff(f, 1, fa, ctx);
	fq_poly_set_coeff(f, 2, f2, ctx);
	fq_poly_set_coeff(f, 3, f3, ctx);
	fq_poly_set_coeff(f, 4, f4, ctx);
	fq_poly_set_coeff(f, 5, f5, ctx);
	fq_poly_set_coeff(f, 6, f6, ctx);
	fq_poly_set_coeff(f, 7, f7, ctx);


	flint_printf("f : ");	
	fq_poly_print(f,ctx);
	flint_printf("\ng : ");
	fq_poly_set_coeff(g, 2, fa, ctx);
	fq_poly_print(g,ctx);
	flint_printf("\nh : ");
	fq_poly_set_coeff(h, 1, fa, ctx);
	fq_poly_set_coeff(h, 8, f1, ctx);
	fq_poly_print(h,ctx);
	flint_printf("\n");

	fq_poly_compose_mod_kedlaya_umans(fg,f,g,h,d,ctx);

	printf("Res\n");

	fq_poly_print(fg,ctx);
	printf("\n");

	printf("Vrai res\n");
	fq_poly_compose_mod(fg2, f, g, h, ctx);

	fq_poly_print(fg2,ctx);
	flint_printf("\n");

	flint_printf("Booléen d'égalité : %d\n",fq_poly_equal(fg,fg2,ctx));


	fmpz_clear(p);
	fmpz_clear(p1);
	fmpz_clear(big);
	fq_poly_clear(f,ctx);
	fq_poly_clear(g,ctx);
	fq_poly_clear(h,ctx);
	fq_poly_clear(fg,ctx);
	fq_poly_clear(fg2,ctx);
	fq_clear(f0,ctx);
	fq_clear(f1,ctx);
	fq_clear(f2,ctx);
	fq_clear(f3,ctx);
	fq_clear(f4,ctx);
	fq_clear(f5,ctx);
	fq_clear(f6,ctx);
	fq_clear(f7,ctx);
	fq_clear(fp,ctx);
	fq_clear(fa,ctx);
	fq_ctx_clear(ctx);

	return 0;
}
