

// Version Fq
void fq_poly_compose_mod_kedlaya_umans(fq_poly_t fg, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, slong d, fq_ctx_t ctx) {

	slong n, i, j, m, dm, N, Nm, degree;
	fmpz* p;
	fmpz_t q;
	fq_multi_poly_t f_prime;
	fq_struct *vect_beta, *vect_alpha, *vect_alpha2, *vect_eval;
	fq_poly_t gi;


	// Vérification des paramètres

	n = FLINT_MAX(fq_poly_length(f,ctx),fq_poly_length(g,ctx));
	i = fq_poly_length(h,ctx);

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

	p = fq_ctx_prime(ctx);
	degree = fq_ctx_degree(ctx);

	m = n_clog(n,d);
	dm = n_pow(d,m);
	N = n_pow(d,m)*m*d;
	Nm = N*m;


	if(fmpz_cmp_ui(p,N) < 0) {
		flint_printf("\nErreur, pas assez de points d'interpolation !\n");
		flint_printf("\np tient sur un mot machine, utilisez donc la version fq_nmod !\n");
		exit(1);
	}

	// Étape 1 //////////////////////////////

	fq_multi_poly_init(f_prime, dm, d, m, ctx);
	_fq_vec_zero(f_prime->poly, dm, ctx);

	for(i = 0 ; i < n ; i++) {
		fq_poly_get_coeff(&(f_prime->poly[i]),f,i,ctx);
	}


	// Étape 3 //////////////////////////////

	vect_beta = _fq_vec_init(N,ctx);

	for(i = 0 ; i < N ; i++) {
		fmpz_poly_set_coeff_ui(vect_beta + i, 0, i);
	}

	vect_alpha = _fq_vec_init(N*m, ctx);


	fq_poly_init(gi,ctx);
	fq_poly_set(gi,g,ctx);
	fq_poly_evaluate_fq_vec(vect_alpha, gi, vect_beta, N, ctx);
	for(i = 1 ; i < m ; i++) {

		// Étape 2 //////////////////////////////
		// Calcul du g^(d^i) avec g^(d^i) = [ (g^(d^(i-1))) ^ d ] mod h
		fq_poly_powmod_ui_binexp(gi,gi,d,h,ctx);

		fq_poly_evaluate_fq_vec(vect_alpha+i*N, gi, vect_beta, N, ctx);
	}

	fq_poly_clear(gi, ctx);


	// Transposition
	// Passe de m lignes N colonnes
	// à N lignes m colonnes ==> nécessaire pour EMM dans Fp

	vect_alpha2 = _fq_vec_init(Nm,ctx);
	for(i = 0 ; i < m ; i++) {
		for(j = 0 ; j < N ; j++) {
			fq_set(vect_alpha2 + j*m + i,vect_alpha + i*N + j,ctx);
		}
	}
	vect_alpha = vect_alpha2;


	// Étape 4 //////////////////////////////
	// EMM

	vect_eval = _fq_vec_init(N,ctx);

	fq_multi_poly_multimodular(vect_eval, N, f_prime, vect_alpha, ctx);

	_fq_vec_clear(vect_alpha,Nm,ctx);
	fq_multi_poly_clear(f_prime,ctx);

	// Étape 5 //////////////////////////////

	fq_poly_interpolate_fq_vec_fast(fg, vect_beta, vect_eval, N, ctx);


	_fq_vec_clear(vect_beta,N,ctx);

	_fq_vec_clear(vect_eval,N,ctx);

	// Étape 6 //////////////////////////////

	fq_poly_rem(fg,fg,h,ctx);

}


