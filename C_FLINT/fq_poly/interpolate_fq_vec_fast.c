

void fq_poly_interpolate_fq_vec_fast(fq_poly_t res, const fq_struct* xs, const fq_struct* ys, slong n, const fq_ctx_t ctx) {

	slong height = FLINT_CLOG2(n);
	slong sup = 1 << height;

	// Construction de l'arbre des sous-produits
    fq_poly_struct ** tree;

    tree = _fq_poly_tree_alloc(n,ctx);
    _fq_poly_tree_build(tree, xs, n, ctx);


	// Calcul de f'
	fq_poly_t f_prime;
	fq_poly_init(f_prime, ctx);

	fq_poly_mul(f_prime,tree[height-1],tree[height-1]+1,ctx);

	fq_poly_derivative(f_prime , f_prime, ctx);


	// Évaluation multipoint en la dérivée
	fq_struct* res2 = _fq_vec_init(n,ctx);

    _fq_poly_evaluate_fq_vec_fast_precomp(res2, f_prime->coeffs, f_prime->length, tree, n, ctx);

	fq_poly_clear(f_prime,ctx);


	fq_struct* res3 = _fq_vec_init(n,ctx);

	slong i, j;
	for(i = 0 ; i < n ; i++) {
		fq_inv(res2+i, res2+i, ctx);
		fq_mul(res3+i, ys+i, res2+i, ctx);
	}

	_fq_vec_clear(res2,n,ctx);



	fq_poly_t* vect_g = malloc(sizeof(fq_poly_t)*sup);
	for(i = 0 ; i < n ; i++) {
		fq_poly_init(vect_g[i],ctx);
		fq_poly_set_fq(vect_g[i],res3+i,ctx);
	}

	_fq_vec_clear(res3,n,ctx);

	for(i = n ; i < sup ; i++) {
		fq_poly_init(vect_g[i],ctx);
		fq_poly_zero(vect_g[i],ctx);
	}


	slong lim = sup;
	slong  nlim = n;
	fq_poly_t var;
	fq_poly_init(var,ctx);
	slong cpt;
	// Parcours de la hauteur de l'arbre
	for(i = 0 ; i < height ; i++) {
		// Parcours de la largeur
		cpt = 0;
		for(j = 0 ; j < lim ; j += 2) {

			if(j+1 >= nlim) {
				fq_poly_set(var,vect_g[j],ctx);
			}
			else {
				fq_poly_mul(var,tree[i]+j+1,vect_g[j],ctx);
			}


			if(j >= nlim) {
				fq_poly_set(vect_g[cpt],vect_g[j+1],ctx);
			}
			else {
				fq_poly_mul(vect_g[cpt],tree[i]+j,vect_g[j+1],ctx);
			}


			fq_poly_add(vect_g[cpt],vect_g[cpt],var,ctx);


			cpt++;
		}
		lim = lim >> 1;
		nlim = (nlim+1) >> 1;
		
	}

    _fq_poly_tree_free(tree, n, ctx);
	fq_poly_clear(var,ctx);

	fq_poly_set(res,vect_g[0],ctx);

	for(i = 0 ; i < sup ; i++) {
		fq_poly_clear(vect_g[i],ctx);
	}

}
