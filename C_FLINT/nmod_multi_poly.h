#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "fmpz.h"


typedef struct
{
    mp_ptr coeffs;
//    slong alloc;
    slong length;
    nmod_t mod;
	// d-1 est le degré max en une variable
	slong d;
	// nombre de variable
	slong m;
} nmod_multi_poly_struct;

typedef nmod_multi_poly_struct nmod_multi_poly_t[1];




void nmod_multi_poly_init(nmod_multi_poly_t poly, slong length, mp_limb_t n, slong d, slong m);
void nmod_multi_poly_clear(nmod_multi_poly_t poly);
