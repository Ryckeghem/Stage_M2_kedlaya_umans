

void fmpz_mod_poly_interpolate_fmpz_vec_fast(fmpz_mod_poly_t res, const fmpz * xs, const fmpz * ys, slong n, const fmpz_t mod) {

	slong i,j, height = FLINT_CLOG2(n), sup = 1 << height, lim = sup, nlim = n, cpt;
	fmpz_mod_poly_t f_prime, cast2, *vect_g, var;
	fmpz_poly_t cast;
	fmpz *res2, *res3;

	// Construction de l'arbre des sous-produits
    fmpz_poly_struct ** tree;

    tree = _fmpz_mod_poly_tree_alloc(n);
    _fmpz_mod_poly_tree_build(tree, xs, n, mod);


	// Calcul de f'
	fmpz_mod_poly_init(f_prime, mod);

	fmpz_poly_init(cast);

	fmpz_poly_mul(cast,tree[height-1],tree[height-1]+1);
	fmpz_mod_poly_set_fmpz_poly(f_prime,cast);

	fmpz_poly_clear(cast);

	fmpz_mod_poly_derivative(f_prime , f_prime);


	// Évaluation multipoint en la dérivée
	res2 = _fmpz_vec_init(n);

    _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(res2, f_prime->coeffs, f_prime->length, tree, n, mod);

	fmpz_mod_poly_clear(f_prime);


	res3 = _fmpz_vec_init(n);

	for(i = 0 ; i < n ; i++) {
		fmpz_invmod(res2+i, res2+i, mod);
		fmpz_mul(res3+i, ys+i, res2+i);
	}

	_fmpz_vec_clear(res2,n);


	fmpz_mod_poly_init(cast2,mod);

	vect_g = malloc(sizeof(fmpz_mod_poly_t)*sup);
	for(i = 0 ; i < n ; i++) {
		fmpz_mod_poly_init(vect_g[i],mod);
		fmpz_mod_poly_set_fmpz(vect_g[i],res3+i);
	}

	_fmpz_vec_clear(res3,n);

	for(i = n ; i < sup ; i++) {
		fmpz_mod_poly_init(vect_g[i],mod);
		fmpz_mod_poly_set_ui(vect_g[i],0);
	}


	fmpz_mod_poly_init(var,mod);
	// Parcours de la hauteur de l'arbre
	for(i = 0 ; i < height ; i++) {
		// Parcours de la largeur
		cpt = 0;
		for(j = 0 ; j < lim ; j += 2) {

			if(j+1 >= nlim) {
				fmpz_mod_poly_set(var,vect_g[j]);
			}
			else {
				fmpz_mod_poly_set_fmpz_poly(cast2,tree[i]+j+1);
				fmpz_mod_poly_mul(var,cast2,vect_g[j]);
			}


			if(j >= nlim) {
				fmpz_mod_poly_set(vect_g[cpt],vect_g[j+1]);
			}
			else {
				fmpz_mod_poly_set_fmpz_poly(cast2,tree[i]+j);
				fmpz_mod_poly_mul(vect_g[cpt],cast2,vect_g[j+1]);
			}

			fmpz_mod_poly_add(vect_g[cpt],vect_g[cpt],var);


			cpt++;
		}
		lim = lim >> 1;
		nlim = (nlim+1) >> 1;
		
	}

    _fmpz_mod_poly_tree_free(tree, n);
	fmpz_mod_poly_clear(cast2);
	fmpz_mod_poly_clear(var);

	fmpz_mod_poly_set(res,vect_g[0]);

	for(i = 0 ; i < sup ; i++) {
		fmpz_mod_poly_clear(vect_g[i]);
	}

}
