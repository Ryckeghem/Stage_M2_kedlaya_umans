

void nmod_multi_poly_fft(mp_ptr vect_points, const mp_ptr f, slong nb_coef, slong p, slong m, mp_limb_t w) {

	slong i, j, pm, nb_coef_m;
	nmod_poly_t poly;


	nmod_poly_init(poly, p);


	if(m == 1) {

		for(i = 0 ; i < nb_coef ; i++) {
			nmod_poly_set_coeff_ui(poly, i, f[i]);
		}

		// On transforme le polynôme en degré p-2, avec X^(p-1) = 1
		if(nb_coef == p) {
			nmod_poly_set_coeff_ui(poly, 0, n_addmod(f[0], f[p-1], p));
			nmod_poly_set_coeff_ui(poly, p-1, 0);
		}

		fft_bluestein2(vect_points, poly, w);

		vect_points[0] = f[0];
	}

	else {

		pm = n_pow(p,m-1);
		nb_coef_m = n_pow(nb_coef,m-1);

		// Récursion sur les coefficients de f
		for(i = 0 ; i < nb_coef ; i++) {
			nmod_multi_poly_fft(vect_points + i*pm, f + i*nb_coef_m, nb_coef, p, m-1, w);
		}

		mp_ptr vect_FFT = _nmod_vec_init(p);

		for(i = 0 ; i < pm ; i++) {
			for(j = 0 ; j < nb_coef ; j++) {
				nmod_poly_set_coeff_ui(poly, j, vect_points[i + j*pm]);
			}


			// On transforme le polynôme en degré p-2, avec X^(p-1) = 1
			if(nb_coef == p) {
				nmod_poly_set_coeff_ui(poly, 0, n_addmod(vect_points[i], vect_points[i + (p-1)*pm], p));
				nmod_poly_set_coeff_ui(poly, p-1, 0);
			}

			fft_bluestein2(vect_FFT, poly, w);
			vect_FFT[0] = vect_points[i];


			for(j = 0 ; j < p ; j++) {
				vect_points[i + j*pm] = vect_FFT[j];
			}

		}

		_nmod_vec_clear(vect_FFT);

	}

	nmod_poly_clear(poly);
}
