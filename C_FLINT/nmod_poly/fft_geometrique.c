

void fft_geometrique(mp_ptr res, const nmod_poly_t g, mp_limb_t w) {

	slong n = nmod_poly_length(g);
	mp_limb_t p = nmod_poly_modulus(g);

	// g(X) = 0
	if(n == 0) {
		_nmod_vec_zero(res,p);
		return;
	}


	res[0] = nmod_poly_get_coeff_ui(g, 0);
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


	mp_limb_t p1 = p-1, p1_demi = p1 >> 1;
	slong i, ti;
	mp_ptr w_vec, t_i;
	nmod_poly_t B, C;


	w_vec = _nmod_vec_init(p1);
	w_vec[0] = 1;
	for(i = 1 ; i < p1 ; i++) {
		w_vec[i] = (w_vec[i-1]*w)%p;
	}


	t_i = _nmod_vec_init(p1);
	t_i[0] = 0;
	nmod_poly_init(B, p);
	nmod_poly_set_coeff_ui(B,0,1);
	nmod_poly_set_coeff_ui(B,p1,w_vec[p1_demi]);
	for(i = 0 ; i < p1-1 ; ) {
		t_i[i+1] = n_addmod(t_i[i],i,p1);
		i++;
		nmod_poly_set_coeff_ui(B,i, w_vec[t_i[i]]);
		nmod_poly_set_coeff_ui(B,i+p1, w_vec[n_addmod(t_i[i],p1_demi,p1)]);
	}


	nmod_poly_init(C, p);
	for(i = 0 ; i < p1 ; i++) {
		res[i+1] = w_vec[n_negmod(t_i[i],p1)];
		nmod_poly_set_coeff_ui(C,p1-i-1, (nmod_poly_get_coeff_ui(g,i)*res[i+1])%p);
	}

	// Mul high
	nmod_poly_mul(B,B,C);
	// Plus lent : nmod_poly_mulhigh(B,B,C,p-2);


	nmod_poly_clear(C);

	for(i = 0 ; i < p1 ; i++) {
		t_i[i] = (res[i+1]*nmod_poly_get_coeff_ui(B,p-2+i))%p;
	}

	nmod_poly_clear(B);

	// On réordonne
	for(i = 0 ; i < p1 ; i++) {
		res[w_vec[i]] = t_i[i];
	}

	_nmod_vec_clear(t_i);
	_nmod_vec_clear(w_vec);
}


