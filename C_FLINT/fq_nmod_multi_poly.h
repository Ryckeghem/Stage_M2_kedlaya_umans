#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "fq.h"
#include "fq_nmod_poly.h"


typedef struct {
	fq_nmod_struct* poly;
//    slong alloc;
    slong length;
	fmpz p;
	// d-1 est le degré max en une variable
	slong d;
	// nombre de variable
	slong m;
} fq_nmod_multi_poly_struct;

typedef fq_nmod_multi_poly_struct fq_nmod_multi_poly_t[1];




void fq_nmod_multi_poly_init(fq_nmod_multi_poly_t f, slong length, slong d, slong m, fq_nmod_ctx_t ctx);
void fq_nmod_multi_poly_clear(fq_nmod_multi_poly_t f, fq_nmod_ctx_t ctx);
