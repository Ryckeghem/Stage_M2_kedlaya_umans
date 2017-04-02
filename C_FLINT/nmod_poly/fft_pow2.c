

// FFT = Évaluation multipoints
// Ne marche qu'à partir d'un point
// p la caractéristique
// w la racine primitive
void nmod_poly_fft_pow2(mp_ptr res, const nmod_poly_t f, mp_limb_t w) {

	slong i, p = nmod_poly_modulus(f);

	slong N = p-1;

	// g(X) = 0
	if(f->length == 0) {
		for(i = 0 ; i < p ; i++) {
			res[i] = 0;
		}
		return;
	}


	// On traite le point 0 qui n'est pas dans le groupe multiplicatif
	res[0] = nmod_poly_get_coeff_ui(f,0);

	if(N == 0) {
		return;
	}

	// Si un seul point
	// Il suffit de sommer les coefficients
	if(N == 1) {
		res[1] = nmod_poly_get_coeff_ui(f,0);
		
		for(i = 1 ; i < f->length; i++) {
			res[1] = n_addmod(res[1], nmod_poly_get_coeff_ui(f,i), 2);
		}
		return;
	}

	// p-1 = 2^k
	mp_limb_t k = n_flog(N,2);

	slong n = N;


	// Arbre d'Évaluation
	slong nmul, nmul2, i2, i3, i4;
	nmod_poly_t* Levaluate = malloc((2*N-1)*sizeof(nmod_poly_t));
	for(i2 = 0 ; i2 < (2*N-1) ; i2++) {
		nmod_poly_init(Levaluate[i2],p);
	}
	nmod_poly_set(Levaluate[0],f);

	nmul = 2;
	nmul2 = 1;
	// i2 parcourt Lpoly
	i2 = 1;
	// i4 parcourt Levaluate
	i4 = 0;

	nmod_poly_t Xmod;
	nmod_poly_init(Xmod, p);

	mp_limb_t w_pow = w, w_pow_i;


	while(nmul <= n) {
		// Pour un niveau de l'arbre
		i3 = 0;


		nmod_poly_set_coeff_ui(Xmod, n/nmul, 1);
		while(i3 <= nmul2-1) {
			// Je calcule le niveau en dessous
			nmod_poly_set_coeff_ui(Xmod, 0, N);
			nmod_poly_rem(Levaluate[i2],Levaluate[i4 + i3],Xmod);
			// Devient F(w*X)
			nmod_poly_set_coeff_ui(Xmod, 0, 1);
			nmod_poly_rem(Levaluate[i2+1],Levaluate[i4 + i3],Xmod);

			w_pow_i = 1;
			for(i = 1 ; i < Levaluate[i2+1]->length ; i++) {
				w_pow_i = (w_pow_i*w_pow)%p;
				// TODO : utiliser nmod_poly_get_coeff_ui et *_set_* serait sans doute plus propre
				Levaluate[i2+1]->coeffs[i] *= w_pow_i;
				Levaluate[i2+1]->coeffs[i] %= p;
				
			}

			i2 += 2;
			i3 += 1;
		}
		nmod_poly_set_coeff_ui(Xmod, n/nmul, 0);
		w_pow = (w_pow*w_pow)%p;
		// Mise à jour du point de départ à réduire modulo ...
		i4 += nmul2;
		nmul2 = nmul;
		nmul = nmul*2;
	}


	// Pour réordonner à la fin

	mp_ptr w_vec = _nmod_vec_init(p-1);
	w_vec[0] = 1;
	for(i = 1 ; i < p-1 ; i++) {
		w_vec[i] = (w_vec[i-1]*w)%p;
	}

	for(i = 1 ; i < p ; i++) {
		if(nmod_poly_length(Levaluate[N -2 + i]) != 0) {
			res[w_vec[n_revbin(i-1,k)]] = nmod_poly_get_coeff_ui(Levaluate[N -2 + i],0);
		}
		else {
			res[w_vec[n_revbin(i-1,k)]] = 0;
		}
	}

	_nmod_vec_clear(w_vec);


}
