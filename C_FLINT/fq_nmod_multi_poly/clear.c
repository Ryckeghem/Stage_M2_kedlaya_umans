

void fq_nmod_multi_poly_clear(fq_nmod_multi_poly_t f, fq_nmod_ctx_t ctx) {
	_fq_nmod_vec_clear(f->poly,f->length,ctx);
}
