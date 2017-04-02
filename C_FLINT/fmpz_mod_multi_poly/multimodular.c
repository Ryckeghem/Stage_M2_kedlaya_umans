

void fmpz_mod_multi_poly_multimodular(fmpz* vect_res, slong N, const fmpz_mod_multi_poly_t multi_f, const fmpz* vect_alpha) {

	slong i, j, k, dm = n_pow(multi_f->d,multi_f->m), Nm = N * multi_f->m;
	fmpz* f = multi_f->poly;
	const mp_limb_t* vect_prime;
	n_primes_t iter;
	fmpz_t var, prod;
	mp_limb_t prime_max;


	fmpz_init(var);

	fmpz_set(var, &(multi_f->p));
	fmpz_sub_ui(var, var, 1);
	fmpz_pow_ui(var, var, (multi_f->m*(multi_f->d -1))+1);
	fmpz_mul_ui(var, var, dm);


	// Calcul de k
	n_primes_init(iter);

	fmpz_init(prod);
	fmpz_set_ui(prod, 1);
	j = 0;

	while(fmpz_cmp(prod, var) < 0) {
		fmpz_mul_ui(prod, prod, n_primes_next(iter));
		j++;
	}
	k = j;

	fmpz_clear(var);
	fmpz_clear(prod);

	n_primes_clear(iter);

	// Vecteur des nombres premiers
	vect_prime = n_primes_arr_readonly(k);

	prime_max = point_fixe(vect_prime[k-1], multi_f->d, multi_f->m, vect_prime);


	mp_ptr f_bar, vect_alpha_prime, vect_eval, vect_tmp;
	nmod_multi_poly_t f_k;
	fmpz_comb_t comb;
	fmpz_comb_temp_t temp;


	// Polynôme qui évalue, version modulée
	f_bar = _nmod_vec_init(dm);
	// Points à évaluer, version modulée
	vect_alpha_prime = _nmod_vec_init(Nm);
	// Vecteur résultat de l'évaluation
	vect_eval = _nmod_vec_init(N*k);


	f_k->length = multi_f->length;
	f_k->d = multi_f->d;
	f_k->m = multi_f->m;	


	for(j = 0 ; j < k ; j++) {


		// Calcul de f_h
		for(i = 0 ; i < dm ; i++) {
			fmpz_mod_ui(var, f+i, vect_prime[j]);
			f_bar[i] = fmpz_get_ui(var);
		}

		f_k->coeffs = f_bar;
		f_k->mod.n = vect_prime[j];

		// Calcul des alpha_h
		for(i = 0 ; i < Nm ; i++) {
			fmpz_mod_ui(var, vect_alpha+i, vect_prime[j]);
			vect_alpha_prime[i] = fmpz_get_ui(var);
		}


		if(vect_prime[j] <= 1.2*prime_max) {
			nmod_multi_poly_evaluate_nmod_vec_fast(vect_eval + j*N, N, f_k, vect_alpha_prime);
		}
		else {
			// Récursion
			nmod_multi_poly_multimodular_primes(vect_eval + j*N, N, f_k, vect_alpha_prime, vect_prime, prime_max);
		}

	}


	fmpz_clear(var);

	// Libère f_bar
	nmod_multi_poly_clear(f_k);
	_nmod_vec_clear(vect_alpha_prime);	


	fmpz_comb_init(comb, vect_prime, k);
	fmpz_comb_temp_init(temp , comb);

	vect_tmp = _nmod_vec_init(k);
	// Transposé
	for(i = 0 ; i < N ; i++) {

		for(j = 0 ; j < k ; j++) {
			vect_tmp[j] = vect_eval[i + j*N];
		}

		fmpz_multi_CRT_ui(vect_res+i, vect_tmp, comb, temp, 0);

		fmpz_mod(vect_res+i,vect_res+i,&(multi_f->p));
	}

	// Le supprimer avant fait planter
	n_cleanup_primes();

	_nmod_vec_clear(vect_eval);
	_nmod_vec_clear(vect_tmp);

	fmpz_comb_clear(comb);
	fmpz_comb_temp_clear(temp);
}




