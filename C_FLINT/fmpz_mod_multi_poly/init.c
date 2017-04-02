
void fmpz_mod_multi_poly_init(fmpz_mod_multi_poly_t f, slong length, const fmpz p, slong d, slong m) {
	f->poly = _fmpz_vec_init(length);
	f->length = length;
	fmpz_set(&(f->p),&p);
	f->d = d;
	f->m = m;
}
