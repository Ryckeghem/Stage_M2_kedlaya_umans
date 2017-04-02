

// On suppose p premier dans fmpz_mod_poly_t ==> corps fini Fp
void fmpz_mod_poly_compose_mod_kedlaya_umans(fmpz_mod_poly_t fg, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_t h, slong d) {

	slong n, i, j, m, dm, N, Nm;
	fmpz_t p;
	fmpz_mod_multi_poly_t f_prime;
	fmpz *vect_beta, *vect_alpha, *vect_alpha2, *vect_eval;
	fmpz_mod_poly_t gi;


	// Vérification des paramètres

	n = FLINT_MAX(fmpz_mod_poly_length(f),fmpz_mod_poly_length(g));
	i = fmpz_mod_poly_length(h);

	if(n >= i) {
		printf("Erreur, f ou g n'a pas été modulé par h.\n");
		exit(1);
	}

	// Les polynômes f et g étant réduits modulo h, n vaut le nombre de coefficients de h
	n = i;

	if(d < 2) {
		printf("d < 2\n");
		exit(1);
	}

	if(d >= n) {
		printf("d >= n\n");
		exit(1);
	}

	fmpz_init(p);
	fmpz_set(p,fmpz_mod_poly_modulus(f));

	m = n_clog(n,d);
	dm = n_pow(d,m);
	N = dm*m*d;
	Nm = N*m;

	// N'arrive pas si le fmpz fait plus de 64 bits
	if(fmpz_cmp_ui(p,N) < 0) {
		flint_printf("\nErreur, pas assez de points d'interpolation !\n");
		exit(1);
	}


	// Étape 1 //////////////////////////////


	fmpz_mod_multi_poly_init(f_prime, dm, *p, d, m);
	_fmpz_vec_zero(f_prime->poly, dm);

	for(i = 0 ; i < n ; i++) {
		fmpz_mod_poly_get_coeff_fmpz(&(f_prime->poly[i]),f,i);
	}


	// Étape 3 //////////////////////////////


	vect_beta = _fmpz_vec_init(N);

	for(i = 0 ; i < N ; i++) {
		fmpz_set_ui(vect_beta + i, i);
	}

	vect_alpha = _fmpz_vec_init(Nm);


	fmpz_mod_poly_init(gi,p);
	fmpz_mod_poly_set(gi,g);
	fmpz_mod_poly_evaluate_fmpz_vec(vect_alpha, gi, vect_beta, N);
	for(i = 1 ; i < m ; i++) {

		// Étape 2 //////////////////////////////

		// Calcul de g^(d^i) avec g^(d^i) = [ (g^(d^(i-1))) ^ d ] mod h
		fmpz_mod_poly_powmod_ui_binexp(gi,gi,d,h);

		fmpz_mod_poly_evaluate_fmpz_vec(vect_alpha+i*N, gi, vect_beta, N);
	}

	fmpz_mod_poly_clear(gi);


	// Transposition
	// Passe de m lignes N colonnes
	// à N lignes m colonnes ==> nécessaire pour EMM dans Fp
	vect_alpha2 = _fmpz_vec_init(Nm);
	for(i = 0 ; i < m ; i++) {
		for(j = 0 ; j < N ; j++) {
			fmpz_set(vect_alpha2 + j*m + i,vect_alpha + i*N + j);
		}
	}
	_fmpz_vec_clear(vect_alpha,Nm);
	vect_alpha = vect_alpha2;


	// Étape 4 //////////////////////////////
	// EMM


	vect_eval = _fmpz_vec_init(N);

	fmpz_mod_multi_poly_multimodular(vect_eval, N, f_prime, vect_alpha);

	_fmpz_vec_clear(vect_alpha,Nm);
	fmpz_mod_multi_poly_clear(f_prime);


	// Étape 5 //////////////////////////////


	fmpz_mod_poly_interpolate_fmpz_vec_fast(fg, vect_beta, vect_eval, N, p);


	_fmpz_vec_clear(vect_beta,N);
	_fmpz_vec_clear(vect_eval,N);
	fmpz_clear(p);


	// Étape 6 //////////////////////////////

	fmpz_mod_poly_rem(fg,fg,h);

}


