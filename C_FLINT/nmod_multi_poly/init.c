

void nmod_multi_poly_init(nmod_multi_poly_t poly, slong length, mp_limb_t n, slong d, slong m) {

	poly->coeffs = (mp_ptr) flint_malloc(length * sizeof(mp_limb_t));;
	poly->length = length;
    poly->mod.n = n;
	poly->d = d;
	poly->m = m;

}
