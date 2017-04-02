
void fq_nmod_multi_poly_init(fq_nmod_multi_poly_t f, slong length, slong d, slong m, fq_nmod_ctx_t ctx) {
	f->poly = _fq_nmod_vec_init(length,ctx);
	f->length = length;
	f->d = d;
	f->m = m;
}
