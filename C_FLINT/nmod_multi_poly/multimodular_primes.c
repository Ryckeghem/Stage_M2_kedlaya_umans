

// Version qui prend la table des nombres premiers en argument
void nmod_multi_poly_multimodular_primes(mp_ptr vect_res, slong N, const nmod_multi_poly_t multi_f, const mp_ptr vect_alpha, const mp_limb_t* vect_prime, mp_limb_t prime_max) {

	slong i, j, k, dm = n_pow(multi_f->d,multi_f->m), Nm = N* multi_f->m;
	mp_ptr f = multi_f->coeffs;
	fmpz_t var, prod;


	fmpz_init(var);

	fmpz_set_ui(var, multi_f->mod.n -1);
	fmpz_pow_ui(var, var, (multi_f->m * (multi_f->d -1))+1);
	fmpz_mul_ui(var, var, dm);


	// Calcul de k
	fmpz_init(prod);
	fmpz_set_ui(prod, 1);
	j = 0;

	while(fmpz_cmp(prod, var) < 0) {
		fmpz_mul_ui(prod, prod, vect_prime[j]);
		j++;
	}
	k = j;

	fmpz_clear(var);
	fmpz_clear(prod);


	// + grand nombre premier <= borne, il n'y a pas de récursion
	if(vect_prime[k-1]  <= 1.2*prime_max) {
		// Cas de base
		nmod_multi_poly_evaluate_nmod_vec_fast(vect_res, N, multi_f, vect_alpha);
		return;
	}

	mp_ptr f_bar, vect_alpha_prime, vect_eval, vect_tmp;
	nmod_multi_poly_t f_k;
	fmpz_comb_t comb;
	fmpz_comb_temp_t temp;
	fmpz_t res_crt;


	// Polynôme qui évalue, version modulée
	f_bar = _nmod_vec_init(dm);
	// Points à évaluer, version modulée
	vect_alpha_prime = _nmod_vec_init(Nm);
	// Vecteur résultat de l'évaluation
	vect_eval = _nmod_vec_init(N*k);

	f_k->length = multi_f->length;
	f_k->d = multi_f->d;
	f_k->m = multi_f->m;	


	// Parcours des nombres premiers
	for(j = 0 ; j < k ; j++) {

		// Calcul de f_h
		for(i = 0 ; i < dm ; i++) {
			f_bar[i] = f[i]%vect_prime[j];
		}

		f_k->coeffs = f_bar;
		f_k->mod.n = vect_prime[j];

		// Calcul des alpha_h
		for(i = 0 ; i < Nm ; i++) {
			vect_alpha_prime[i] = vect_alpha[i]%vect_prime[j];
		}


		if(vect_prime[j] <= 1.2*prime_max) {
			nmod_multi_poly_evaluate_nmod_vec_fast(vect_eval + j*N, N, f_k, vect_alpha_prime);
		}
		else {
			// Récursion
			nmod_multi_poly_multimodular_primes(vect_eval + j*N, N, f_k, vect_alpha_prime, vect_prime, prime_max);
		}

	}

	// Libère f_bar
	nmod_multi_poly_clear(f_k);
	_nmod_vec_clear(vect_alpha_prime);	

	fmpz_comb_init(comb, vect_prime, k);
	fmpz_comb_temp_init(temp , comb);

	vect_tmp = _nmod_vec_init(k);

	fmpz_init(res_crt);


	for(i = 0 ; i < N ; i++) {

		for(j = 0 ; j < k ; j++) {
			vect_tmp[j] = vect_eval[i + j*N];
		}

		fmpz_multi_CRT_ui(res_crt, vect_tmp, comb, temp, 0);
		fmpz_mod_ui(res_crt,res_crt,multi_f->mod.n);
		vect_res[i] = fmpz_get_ui(res_crt);
	}

	// Le supprimer avant fait planter
	_nmod_vec_clear(vect_eval);
	_nmod_vec_clear(vect_tmp);


	fmpz_clear(res_crt);
	fmpz_comb_clear(comb);
	fmpz_comb_temp_clear(temp);

}

