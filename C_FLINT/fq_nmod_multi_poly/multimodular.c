

// N'est pas valide dans certains cas, c'est probablement fq_nmod_multi_poly_multimodular qui ne marche pas bien pour les extensions (problème de cast de fq_nmod à fq)
void fq_nmod_multi_poly_multimodular(fq_nmod_struct* vect_res, slong N, const fq_nmod_multi_poly_t multi_f, const fq_nmod_struct* vect_alpha, fq_nmod_ctx_t ctx) {

	slong i,j,l,Nm;
	fq_struct *vect_res2, *vect_alpha2;
	fq_multi_poly_t multi_f2;
	// Création d'un contexte cast de fq
	fq_ctx_t ctx2;
	fmpz_mod_poly_t mod;
	fmpz* p = fq_nmod_ctx_prime(ctx);


	fmpz_mod_poly_init(mod, p);
	l = nmod_poly_length(ctx->modulus);
	for(j = 0 ; j < l ; j++) {
		fmpz_mod_poly_set_coeff_ui(mod, j, (ctx->modulus)->coeffs[j]);
	}

	fq_ctx_init_modulus(ctx2, mod, "b");
	fmpz_mod_poly_clear(mod);


	fq_multi_poly_init(multi_f2, n_pow(multi_f->d,multi_f->m), multi_f->d, multi_f->m, ctx2);
	for(i = 0 ; i < multi_f->length ; i++) {

		l = nmod_poly_length(multi_f->poly+i);
		for(j = 0 ; j < l ; j++) {
			fmpz_poly_set_coeff_ui((multi_f2->poly+i), j, (multi_f->poly+i)->coeffs[j]);
		}
	}

	printf("\nfbarrr :\n");
	_fq_nmod_vec_print(multi_f->poly, multi_f->length, ctx);
	printf("\n");

	_fq_vec_print(multi_f2->poly, multi_f2->length, ctx2);
	printf("\n");



	Nm = N*(multi_f->m);
	vect_alpha2 = _fq_vec_init(Nm,ctx2);
	for(i = 0 ; i < Nm ; i++) {
		l = nmod_poly_length(vect_alpha+i);
		for(j = 0 ; j < l ; j++) {
			fmpz_poly_set_coeff_ui(vect_alpha2+i, j, (vect_alpha+i)->coeffs[j]);
		}
	}


	printf("\nalphaa :\n");
	_fq_nmod_vec_print(vect_alpha, Nm, ctx);
	printf("\n");

	printf("\nalphaa2 :\n");
	_fq_vec_print(vect_alpha2, Nm, ctx2);
	printf("\n");


	vect_res2 = _fq_vec_init(N,ctx2);
	fq_multi_poly_multimodular(vect_res2, N, multi_f2, vect_alpha2, ctx2);

	_fq_vec_clear(vect_alpha2,Nm,ctx2);
	fq_multi_poly_clear(multi_f2, ctx2);


	for(i = 0 ; i < N ; i++) {
		l = fmpz_poly_length(vect_res2+i);
		for(j = 0 ; j < l ; j++) {
			nmod_poly_set_coeff_ui(vect_res+i, j, fmpz_get_ui(&(vect_res2+i)->coeffs[j]));
		}
	}

	_fq_vec_clear(vect_res2,N,ctx2);
	fq_ctx_clear(ctx2);

}


