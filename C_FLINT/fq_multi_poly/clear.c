

void fq_multi_poly_clear(fq_multi_poly_t f, fq_ctx_t ctx) {
	_fq_vec_clear(f->poly,f->length,ctx);
}
