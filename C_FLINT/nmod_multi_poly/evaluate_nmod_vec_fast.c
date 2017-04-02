

void nmod_multi_poly_evaluate_nmod_vec_fast(mp_ptr vect_res, slong N, const nmod_multi_poly_t multi_f, const mp_ptr vect_alpha) {

	slong i,j,pi,nb_coef,*P;
	mp_limb_t w;
	mp_ptr f = multi_f->coeffs, f_bar, vect_points, vect_indice;


	P = (slong*)malloc(sizeof(slong)* multi_f->m);
	P[0] = 1;
	for(pi = 1 ; pi < multi_f->m ; pi++) {
		P[pi] = P[pi-1]* multi_f->mod.n;
	}

	// Modulo X_i ^p - X
	if(multi_f->mod.n < multi_f->d) {

		nb_coef = multi_f->mod.n;
		f_bar = _nmod_vec_init(n_pow(nb_coef,multi_f->m));
		nmod_multi_poly_mod_fermat(f_bar, f, multi_f->mod.n, multi_f->d, multi_f->m, P);

	}
	else {
		nb_coef = multi_f->d;
		f_bar = _nmod_vec_init(n_pow(nb_coef,multi_f->m));
		for(i = 0 ; i < n_pow(multi_f->d,multi_f->m) ; i++) {
			f_bar[i] = f[i];
		}
	}

	vect_points = _nmod_vec_init(n_pow(multi_f->mod.n,multi_f->m));

	w = n_primitive_root_prime(multi_f->mod.n);
	nmod_multi_poly_fft(vect_points, f_bar, nb_coef, multi_f->mod.n, multi_f->m, w);

	_nmod_vec_clear(f_bar);



	vect_indice = _nmod_vec_init(multi_f->m);

	for(i = 0 ; i < N ; i++) {

		// Recherche du point sous forme (Fp)^m
		for(j = 0 ; j < multi_f->m ; j++) {
			vect_indice[j] = vect_alpha[i * multi_f->m + j];
		}

		vect_res[i] = vect_points[Coord(P, vect_indice, multi_f->m)];

	}


	free(P);
	_nmod_vec_clear(vect_indice);
	_nmod_vec_clear(vect_points);

}


