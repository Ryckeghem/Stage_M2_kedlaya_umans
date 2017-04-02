

// Fq_nmod
// N'est pas valide dans certains cas, c'est probablement fq_nmod_multi_poly_multimodular qui ne marche pas bien pour les extensions (problème de cast de fq_nmod à fq)
void fq_nmod_poly_compose_mod_kedlaya_umans(fq_nmod_poly_t fg, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, slong d, fq_nmod_ctx_t ctx) {

	slong n;

	n = FLINT_MAX(fq_nmod_poly_length(f,ctx),fq_nmod_poly_length(g,ctx));
	n = FLINT_MAX(n,fq_nmod_poly_length(h,ctx));

	if(d < 2) {
		printf("d < 2\n");
		exit(1);
	}

	if(d >= n) {
		printf("d >= n\n");
		exit(1);
	}

	slong i;
	slong m, N;
	slong degree = fq_nmod_ctx_degree(ctx);
	fmpz* p;

	p = fq_nmod_ctx_prime(ctx);

	mp_limb_t pp = fmpz_get_ui(p);

	m = n_clog(n,d);

	fmpz_t q;
	fmpz_init(q);

	fmpz_pow_ui(q,p,degree);

	// Étape 1 //////////////////////////////
	slong dm = n_pow(d,m);

	fq_nmod_multi_poly_t f_prime;
	fq_nmod_multi_poly_init(f_prime, dm, d, m, ctx);
	_fq_nmod_vec_zero(f_prime->poly, dm, ctx);

	for(i = 0 ; i < n ; i++) {
		fq_nmod_poly_get_coeff(&(f_prime->poly[i]),f,i,ctx);
	}



	// Étape 3 //////////////////////////////

	N = n_pow(d,m)*m*d;
	fq_nmod_struct* vect_beta = _fq_nmod_vec_init(N,ctx);

	flint_printf("\n");
	flint_printf("n %d d %d m %d N %d p^e ",n,d,m,N);
	fmpz_print(p);
	flint_printf("^%d\n",degree);


	if(fmpz_cmp_ui(q,N) < 0) {
		flint_printf("\nErreur, pas assez de points d'interpolation !\n");
		exit(1);
	}


	mp_ptr cpt_vec = _nmod_vec_init(degree);
	_nmod_vec_zero(cpt_vec,degree);

	slong j;
	if(N < pp) {
		for(i = 0 ; i < N ; i++) {
			nmod_poly_set_coeff_ui(vect_beta + i, 0, i);
		}
	}
	else{

		// Initialisation du premier point à 0
		nmod_poly_set_coeff_ui(vect_beta, 0, 0);

		for(i = 1 ; i < N ; i++) {
			cpt_vec[0] += 1;		

			// si j'ai une retenue
			if(cpt_vec[0] == pp) {

				cpt_vec[0] = 0;
				cpt_vec[1] += 1;

				j = 1;
				while(cpt_vec[j] == pp) {
					cpt_vec[j+1] += 1;
					cpt_vec[j] = 0;
					j++;
				}
			}

			for(j = 0 ; j < degree ; j++) {
				nmod_poly_set_coeff_ui(vect_beta+i, j, cpt_vec[j]);
			}
	//		_fmpz_vec_print(cpt_vec,degree);
	//		flint_printf("\n");
		}
	}
	_nmod_vec_clear(cpt_vec);

	_fq_nmod_vec_print(vect_beta,N,ctx);


	fq_nmod_struct* vect_alpha = _fq_nmod_vec_init(N*m, ctx);


	fq_nmod_poly_t gi;
	fq_nmod_poly_init(gi,ctx);
	fq_nmod_poly_set(gi,g,ctx);
	fq_nmod_poly_evaluate_fq_nmod_vec(vect_alpha, gi, vect_beta, N, ctx);

	flint_printf("\n");
	_fq_nmod_vec_print(vect_alpha,N,ctx);
	flint_printf("\n");

	for(i = 1 ; i < m ; i++) {

		// Étape 2 //////////////////////////////
		// Calcul du g^(d^i) avec g^(d^i) = [ (g^(d^(i-1))) ^ d ] mod h
		fq_nmod_poly_powmod_ui_binexp(gi,gi,d,h,ctx);

		flint_printf("\n");
		fq_nmod_poly_print(gi,ctx);
		flint_printf("\n");

		fq_nmod_poly_evaluate_fq_nmod_vec(vect_alpha+i*N, gi, vect_beta, N, ctx);

		flint_printf("\n");
		_fq_nmod_vec_print(vect_alpha+i*N,N,ctx);
		flint_printf("\n");
	}

	fq_nmod_poly_clear(gi, ctx);

	flint_printf("\n");

	_fq_nmod_vec_print(vect_alpha,N*m,ctx);
	flint_printf("\n");

	// Transposition
	// Passe de m lignes N colonnes
	// à N lignes m colonnes ==> nécessaire pour EMM dans Fp

	fq_nmod_struct* vect_alpha2 = _fq_nmod_vec_init(N*m,ctx);
	for(i = 0 ; i < m ; i++) {
		for(j = 0 ; j < N ; j++) {
			fq_nmod_set(vect_alpha2 + j*m + i,vect_alpha + i*N + j,ctx);
		}
	}
	vect_alpha = vect_alpha2;

	_fq_nmod_vec_print(vect_alpha,N*m,ctx);

	// Étape 4 //////////////////////////////
	// EMM

	flint_printf("\nEtape 4 \n");

	fq_nmod_struct* vect_eval = _fq_nmod_vec_init(N,ctx);

	fq_nmod_multi_poly_multimodular(vect_eval, N, f_prime, vect_alpha, ctx);

	// TODO : Nm = N*m
	_fq_nmod_vec_clear(vect_alpha,N*m,ctx);
	fq_nmod_multi_poly_clear(f_prime,ctx);

	// Étape 5 //////////////////////////////

	flint_printf("\nEtape 5 \n");

	fq_nmod_poly_interpolate_fq_nmod_vec_fast(fg, vect_beta, vect_eval, N, ctx);


	_fq_nmod_vec_clear(vect_beta,N,ctx);

	_fq_nmod_vec_print(vect_eval,N,ctx);
	flint_printf("\n");

	// Étape 6 //////////////////////////////

	flint_printf("\nEtape 6 \n");

	fq_nmod_poly_rem(fg,fg,h,ctx);

}


