

void fmpz_mod_multi_poly_clear(fmpz_mod_multi_poly_t f) {
	_fmpz_vec_clear(f->poly,f->length);
}
