

// On suppose p premier dans fmpz_mod_poly_t ==> corps fini Fp
void nmod_poly_compose_mod_kedlaya_umans(nmod_poly_t fg, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h, slong d) {

	slong n, i, j, m, dm, N, Nm;
	mp_limb_t p;
	nmod_multi_poly_t f_prime;
	mp_ptr vect_beta, vect_alpha, vect_alpha2, vect_eval;
	nmod_poly_t gi;


	// Vérification des paramètres

	n = FLINT_MAX(nmod_poly_length(f),nmod_poly_length(g));
	i = nmod_poly_length(h);

	if(n >= i) {
		printf("Erreur, f ou g n'a pas été modulé par h.\n");
		exit(1);
	}

	// Les polynômes f et g étant modulés par h, n vaut le nombre de coefficients de h
	n = i;

	if(d < 2) {
		printf("d < 2\n");
		exit(1);
	}

	if(d >= n) {
		printf("d >= n\n");
		exit(1);
	}

	p = nmod_poly_modulus(f);

	m = n_clog(n,d);
	dm = n_pow(d,m);
	N = dm*m*d;
	Nm = N*m;


	if(p < N) {
		flint_printf("\nErreur, pas assez de points d'interpolation !\n");
		exit(1);
	}



	// Étape 1 //////////////////////////////
	nmod_multi_poly_init(f_prime, dm, p, d, m);
	_nmod_vec_zero(f_prime->coeffs, dm);

	for(i = 0 ; i < n ; i++) {
		f_prime->coeffs[i] = nmod_poly_get_coeff_ui(f,i);
	}

	// Étape 3 //////////////////////////////

	vect_beta = _nmod_vec_init(N);

	for(i = 0 ; i < N ;) {
		vect_beta[i] = i++;
	}

	vect_alpha = _nmod_vec_init(N*m);


	nmod_poly_init(gi,p);
	nmod_poly_set(gi,g);
	nmod_poly_evaluate_nmod_vec(vect_alpha, gi, vect_beta, N);

	for(i = 1 ; i < m ; i++) {

		// Étape 2 //////////////////////////////
		// Calcul du g^(d^i) avec g^(d^i) = [ (g^(d^(i-1))) ^ d ] mod h
		nmod_poly_powmod_ui_binexp(gi,gi,d,h);

		nmod_poly_evaluate_nmod_vec(vect_alpha+i*N, gi, vect_beta, N);
	}

	nmod_poly_clear(gi);


	// Transposition
	// Passe de m lignes N colonnes
	// à N lignes m colonnes ==> nécessaire pour EMM dans Fp
	vect_alpha2 = _nmod_vec_init(Nm);
	for(i = 0 ; i < m ; i++) {
		for(j = 0 ; j < N ; j++) {
			vect_alpha2[j*m + i] = vect_alpha[i*N + j];
		}
	}
	_nmod_vec_clear(vect_alpha);
	vect_alpha = vect_alpha2;


	// Étape 4 //////////////////////////////
	// EMM

	vect_eval = _nmod_vec_init(N);

	nmod_multi_poly_multimodular(vect_eval, N, f_prime, vect_alpha);

	_nmod_vec_clear(vect_alpha);
	nmod_multi_poly_clear(f_prime);

	// Étape 5 //////////////////////////////

	nmod_poly_interpolate_nmod_vec_fast(fg, vect_beta, vect_eval, N);


	_nmod_vec_clear(vect_beta);
	_nmod_vec_clear(vect_eval);


	// Étape 6 //////////////////////////////

	nmod_poly_rem(fg,fg,h);

}


