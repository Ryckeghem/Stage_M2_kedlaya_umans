

// Considérer n <= p-1, c'est-à-dire g de degré p-2 au plus
void fft_bluestein(mp_ptr res, const nmod_poly_t g, mp_limb_t w) {

	slong i, n, n_low;
	mp_limb_t p = nmod_poly_modulus(g);

	n = nmod_poly_length(g);

	// g(X) = 0
	if(n == 0) {
		for(i = 0 ; i < p ; i++) {
			res[i] = 0;
		}
		return;
	}


	res[0] = nmod_poly_get_coeff_ui(g, 0);
	if(p == 1) {
		return;
	}


	res[1] = 0;
	for(i = 0 ; i < n ; i++) {
		res[1] = n_addmod(res[1],nmod_poly_get_coeff_ui(g, i),p);
	}

	if(p == 2) {
		return;
	}

	if(n > p) {
		flint_printf("Erreur, il faut appliquer X^p = X.\n");
		exit(1);
	}

	if(n == p) {
		flint_printf("Erreur, il faut appliquer X^(p-1) = 1.\n");
		exit(1);
	}


	fq_nmod_ctx_t ctx;
	nmod_poly_t modulus;
	slong i2;

	mp_ptr w_vec;
	fq_nmod_struct* w_fq_vec, *wq_inv;
	fq_nmod_t mul;
	fq_nmod_poly_t f, xq, wq, res_low, res_high;



	// Création d'une extension, on ajoute une racine carré de w
	// X^2 - w = 0
	nmod_poly_init(modulus,p);
	nmod_poly_set_coeff_ui(modulus,0,p-w);
	nmod_poly_set_coeff_ui(modulus,2,1);
	fq_nmod_ctx_init_modulus(ctx , modulus, "a");


	// Cast de g dans l'extension = f
	fq_nmod_poly_init(f,ctx);
	fq_nmod_t cast;
	fq_nmod_init(cast, ctx);

	for(i = 0 ; i < n ; i++) {
		fq_nmod_set_ui(cast, nmod_poly_get_coeff_ui(g, i), ctx);
		fq_nmod_poly_set_coeff(f, i, cast, ctx);
	}



	w_vec = _nmod_vec_init(p-1);
	w_vec[0] = 1;
	for(i = 1 ; i < p-1 ; i++) {
		w_vec[i] = (w_vec[i-1]*w)%p;
	}

	// Tableau cast de w_vec
	w_fq_vec = _fq_nmod_vec_init(p-1,ctx);
	for(i = 0 ; i < p-1 ; i++) {
		nmod_poly_set_coeff_ui(w_fq_vec+i, 0, w_vec[i]);
	}

	fq_nmod_init(mul,ctx);
	fq_nmod_poly_init(xq,ctx);
	fq_nmod_poly_init(wq,ctx);

	wq_inv = _fq_nmod_vec_init(p-1,ctx);


	// Astuce pour avoir la taille max du poly, et éviter le bug du coef
	// qui ne se met pas à jour s'il n'existe pas (0 = n'existe pas)
	// Initialisé à 1 (= w_fq_vec), faut remettre 0 ensuite

	// On fait exprès de dépasser de 1
	fq_nmod_poly_set_coeff(wq,p-1,w_fq_vec,ctx);
	// n-1 <= p-2
	fq_nmod_poly_set_coeff(xq,n,w_fq_vec,ctx);


	fq_nmod_poly_set_coeff(wq,0,w_fq_vec,ctx);

	fq_nmod_set_ui(wq_inv, 1, ctx);

	// Je sais que wq_inv = 1, donc 1*f->coeffs = f->coeffs
	fq_nmod_poly_set_coeff(xq,0,f->coeffs,ctx);


	// Cas pair, on reste dans le corps de base
	for(i = 2 ; i < p-1 ; i += 2) {
		i2 = (i*i >> 1)%(p-1);
		fq_nmod_poly_set_coeff(wq,i,w_fq_vec+i2,ctx);

		// On sait que le nombre est négatif, donc suffit d'ajouter p-1
		nmod_poly_set_coeff_ui(wq_inv+i, 0, w_vec[(p-1-i2)%(p-1)]);
	}


	for(i = 2 ; i < n ; i += 2) {
		fq_nmod_mul(mul, wq_inv+i, f->coeffs+i, ctx);
		fq_nmod_poly_set_coeff(xq,i,mul,ctx);
	}


	_fq_nmod_vec_clear(w_fq_vec,p-1,ctx);

	// Cas impair, on est "imaginaire pur"
	for(i = 1 ; i < p-1 ; i += 2) {
		i2 = (i*i >> 1)%(p-1);
		nmod_poly_set_coeff_ui(wq->coeffs+i,1,w_vec[i2]);

		// - 1 + p-1 = p-2, le -1 est pour passer de puissance -1/2 à puissance 1.2
		nmod_poly_set_coeff_ui(wq_inv+i, 1, w_vec[(p-2-i2)]);

	}

	for(i = 1 ; i < n ; i += 2) {
		fq_nmod_mul(mul, wq_inv+i, f->coeffs+i, ctx);
		fq_nmod_poly_set_coeff(xq,i,mul,ctx);
	}


	fq_nmod_clear(mul,ctx);

	// On écrase le surplus ajouté
	fq_nmod_set_ui(wq->coeffs+p-1,0,ctx);
	fq_nmod_set_ui(xq->coeffs+n,0,ctx);


	fq_nmod_poly_init(res_low,ctx);
	fq_nmod_poly_init(res_high,ctx);

	fq_nmod_poly_mullow(res_low,xq,wq,p-1,ctx);

	fq_nmod_poly_reverse(xq,xq,p-1,ctx);
	fq_nmod_zero(wq->coeffs,ctx);


	fq_nmod_poly_mullow(res_high,xq,wq,p-2,ctx);
	fq_nmod_poly_clear(xq,ctx);
	fq_nmod_poly_clear(wq,ctx);

	fq_nmod_poly_reverse(res_high,res_high,p-1,ctx);

	fq_nmod_poly_add(res_low,res_low,res_high,ctx);
	fq_nmod_poly_clear(res_high,ctx);

	n_low = FLINT_MIN(fq_nmod_poly_length(res_low,ctx),p-1);

	for(i = 0 ; i < n_low ; i++) {
		fq_nmod_mul(wq_inv+i, wq_inv+i, res_low->coeffs+i, ctx);
	}
	fq_nmod_poly_clear(res_low,ctx);


	// On réordonne
	for(i = 1 ; i < p-1 ; i++) {

		if(nmod_poly_length(wq_inv+i) != 0) {
			res[w_vec[p-1-i]] = (wq_inv+i)->coeffs[0];
		}
		else {
			res[w_vec[p-1-i]] = 0;
		}
	}

	_nmod_vec_clear(w_vec);
	_fq_nmod_vec_clear(wq_inv,p-1,ctx);

}
