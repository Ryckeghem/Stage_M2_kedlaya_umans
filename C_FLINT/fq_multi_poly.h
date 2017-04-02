#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "fq.h"
#include "fq_nmod_poly.h"


typedef struct {
	fq_struct* poly;
//    slong alloc;
    slong length;
	fmpz p;
	// d-1 est le degré max en une variable
	slong d;
	// nombre de variable
	slong m;
} fq_multi_poly_struct;

typedef fq_multi_poly_struct fq_multi_poly_t[1];




void fq_multi_poly_init(fq_multi_poly_t f, slong length, slong d, slong m, fq_ctx_t ctx);
void fq_multi_poly_clear(fq_multi_poly_t f, fq_ctx_t ctx);
