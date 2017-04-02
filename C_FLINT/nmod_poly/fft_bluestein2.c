

// Considérer n <= p-1, c'est-à-dire g de degré p-2 au plus
void fft_bluestein2(mp_ptr res, const nmod_poly_t g, mp_limb_t w) {

	slong n = nmod_poly_length(g);
	mp_limb_t p = nmod_poly_modulus(g);

	// g(X) = 0
	if(n == 0) {
		_nmod_vec_zero(res,p);
		return;
	}


	res[0] = nmod_poly_get_coeff_ui(g, 0);
	// Peut-être impossible sous FLINT, res de taille 1
	if(p == 1) {
		return;
	}

	// g(X) = constante
	// Traite le cas p = 2
	if(n == 1) {
		slong i;
		for(i = 1 ; i < p ; i++) {
			res[i] = res[0];
		}
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

	if(p < 20) {
		slong i;
		mp_ptr vect = _nmod_vec_init(p);

		for(i = 1 ; i < p ; i++) {
			vect[i] = i;
		}

		res[0] = g->coeffs[0];
		nmod_poly_evaluate_nmod_vec(res+1,g,vect+1,p-1);

		_nmod_vec_clear(vect);
		return;
	}

	mp_limb_t p1 = p-1, p1_demi = p1 >> 1;
	slong i, i_mod, var, pow_theta = p1_demi%2;
	mp_ptr w_vec, fi, fi_prime;
	nmod_poly_t F, G, FG;



	w_vec = _nmod_vec_init(p1);
	w_vec[0] = 1;
	for(i = 1 ; i < p1 ; i++) {
		w_vec[i] = (w_vec[i-1]*w)%p;
	}


	// La taille de fi devrait être p1_demi+1, mais on garde plus de place
	// pour stocker le résultat final ordonné selon w^i
	fi = _nmod_vec_init(p1);
	fi_prime = _nmod_vec_init(p1_demi);

	i_mod = 0;
	fi[0] = 1;
	fi_prime[0] = 1;

	nmod_poly_init(G, p);
	nmod_poly_set_coeff_ui(G,0,2);
	for(i = 1 ; i < p1_demi ; i++) {

		// i^2 mod p-1
		i_mod = i_mod+(i << 1)-1;
		if(i_mod >= p1) {
			i_mod -= p1;
		}
		fi[i] = w_vec[i_mod];

		
		// i^2 + i mod p-1
		n = i_mod + i;
		if(n >= p1) {
			n -= p1;
		}
		fi_prime[i] = w_vec[n];


		// w^(-i^2) + w^(-i^2-i))
		if(n != 0) {
			n = p1-n;
		}

		var = -i_mod;
		if(var != 0) {
			var += p1;
		}

		nmod_poly_set_coeff_ui(G,i,n_addmod(w_vec[var], w_vec[n], p));

		// Utilisation de la "périodicité" des g_i
		// theta*(w^(-i^2) - w^(-i^2-i))
		if(pow_theta == 0) {
			// theta = 1
			nmod_poly_set_coeff_ui(G,i+p1_demi,n_submod(w_vec[var], w_vec[n], p));
		}
		else {
			// theta = -1
			nmod_poly_set_coeff_ui(G,i+p1_demi,n_submod(w_vec[n], w_vec[var], p));
		}
	}

	// Cas i = p1_demi
	i_mod = i_mod+(i << 1)-1;
	if(i_mod >= p1) {
		i_mod -= p1;
	}
	fi[i] = w_vec[i_mod];


	nmod_poly_set_coeff_ui(G,i,0);


	nmod_poly_init(F, p);

	for(i = 0 ; i <= p1_demi ; i++) {
		nmod_poly_set_coeff_ui(F,i,fi[i]*nmod_poly_get_coeff_ui(g,i));
	}

	for(i = 1 ; i < p1_demi ; i++) {
		// On utilise la symétrie de i^2 mod p-1
		nmod_poly_set_coeff_ui(F,p1_demi+i,fi[p1_demi-i]*nmod_poly_get_coeff_ui(g,p1_demi+i));
	}


	nmod_poly_init(FG, p);

	nmod_poly_mul(FG,F,G);

	nmod_poly_shift_right(G,FG,p1);
	nmod_poly_truncate(FG,p1);
	nmod_poly_add(F, FG, G);


	nmod_poly_scalar_mul_nmod(F, F, (p+1) >> 1);

	nmod_poly_shift_right(G,F,p1_demi);
	nmod_poly_truncate(F,p1_demi);


	nmod_poly_add(FG, F, G);
	nmod_poly_sub(G, F, G);

	nmod_poly_clear(F);

	// theta
	if(pow_theta == 0) {
		// theta = 1
		// La décrémentation permet de ne pas avoir de conflit sur fi
		for(i = p1_demi-1 ; i >= 0 ; i--) {
			fi[i << 1] = (nmod_poly_get_coeff_ui(FG, i)*fi[i])%p;
		}

		for(i = 0 ; i < p1_demi ; i++) {
			fi[(i << 1)+1] = (nmod_poly_get_coeff_ui(G, i)*fi_prime[i])%p;
		}

	}
	else {
		// theta = -1
		// La décrémentation permet de ne pas avoir de conflit sur fi
		for(i = p1_demi-1 ; i >= 0 ; i--) {
			fi[i << 1] = (nmod_poly_get_coeff_ui(G, i)*fi[i])%p;
		}

		for(i = 0 ; i < p1_demi ; i++) {
			fi[(i << 1)+1] = (nmod_poly_get_coeff_ui(FG, i)*fi_prime[i])%p;
		}
	}

	nmod_poly_clear(G);
	nmod_poly_clear(FG);
	_nmod_vec_clear(fi_prime);

	// On réordonne

	for(i = 0 ; i < p1 ; i++) {
		res[w_vec[i]] = fi[i];
	}

	_nmod_vec_clear(w_vec);
	_nmod_vec_clear(fi);

}

