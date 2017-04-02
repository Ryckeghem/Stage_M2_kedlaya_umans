#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"


typedef struct {
	fmpz* poly;
//    slong alloc;
    slong length;
	fmpz p;
	// d-1 est le degré max en une variable
	slong d;
	// nombre de variable
	slong m;
} fmpz_mod_multi_poly_struct;

typedef fmpz_mod_multi_poly_struct fmpz_mod_multi_poly_t[1];




void fmpz_mod_multi_poly_init(fmpz_mod_multi_poly_t f, slong length, fmpz p, slong d, slong m);
void fmpz_mod_multi_poly_clear(fmpz_mod_multi_poly_t f);
