

void fq_multi_poly_multimodular(fq_struct* vect_res, slong N, const fq_multi_poly_t multi_f, const fq_struct* vect_alpha, fq_ctx_t ctx) {

	slong i, j, k, Nm;

	slong e = fq_ctx_degree(ctx);
	fmpz* p = fq_ctx_prime(ctx);
	fmpz_t M, r_prime;


	fmpz_init(M);

	fmpz_set(M, p);
	fmpz_sub_ui(M, M, 1);
	fmpz_mul_ui(M, M, e);
	fmpz_pow_ui(M, M, ((multi_f->d)-1)*(multi_f->m) +1);
	fmpz_mul_ui(M, M, n_pow(multi_f->d,multi_f->m));	
	fmpz_add_ui(M, M, 1);


	k = (e-1)*(multi_f->d)*(multi_f->m);
	fmpz_init(r_prime);

	fmpz_pow_ui(r_prime, M, k+1);


	fmpz_mod_multi_poly_t f_bar;
	fmpz_mod_multi_poly_init(f_bar, multi_f->length, *r_prime, multi_f->d, multi_f->m);

	// Correct
	for(i = 0 ; i < multi_f->length ; i++) {

		// Sustitution de Kronecker
		// Évaluation d'un élément de Fq_mod = fmpz_poly_t en M
		fmpz_poly_evaluate_fmpz(f_bar->poly+i, multi_f->poly+i, M);

	}

	// Points à évaluer, version Kronecker
	Nm = N*(multi_f->m);
	fmpz* vect_alpha_bar = _fmpz_vec_init(Nm);


	fmpz_t var;
	fmpz_init(var);


	for(i = 0 ; i < Nm ; i++) {

		// Sustitution de Kronecker
		// Évaluation d'un élément de Fq = fmpz_poly_t en M
		fmpz_poly_evaluate_fmpz(vect_alpha_bar+i, vect_alpha+i, M);

	}


	// Vecteur résultat de l'évaluation
	fmpz* vect_eval = _fmpz_vec_init(N);

	fmpz_mod_multi_poly_multimodular(vect_eval, N, f_bar, vect_alpha_bar);

	fmpz_mod_multi_poly_clear(f_bar);
	_fmpz_vec_clear(vect_alpha_bar, N);


	fmpz_t q, r;
	fmpz_init(q);
	fmpz_init(r);


	for(i = 0 ; i < N ; i++) {

		fmpz_set(q, vect_eval+i);

		j = 0;
		while(fmpz_is_zero(q) == 0) {
			fmpz_fdiv_qr(q, r, q, M);
			fmpz_mod(r,r,p);
			fmpz_poly_set_coeff_fmpz(vect_res+i,j,r);

			j++;
		}

		while(j < k) {
			fmpz_poly_set_coeff_ui(vect_res+i,j,0);
			j++;
		}

	}

	fmpz_clear(q);
	fmpz_clear(r);
	_fmpz_vec_clear(vect_eval, N);

	// modulo mod r et E(z)
	for(i = 0 ; i < N ; i++) {
		fq_reduce(vect_res+i, ctx);
	}

}
