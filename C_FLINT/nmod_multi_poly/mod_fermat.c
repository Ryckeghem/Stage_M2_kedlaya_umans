
slong Coord(const slong* D, const slong* pos, slong m) {

	slong i, res = 0;

	for(i = 0 ; i < m ; i++) {
		res += pos[i]*D[i];
	}

	return res;
}


slong Coord_mod(const slong* P, slong p, const slong* pos, slong m) {

	slong i, res = 0;
	slong pos_tmp;

	if(p == 2) {
		for(i = 0 ; i < m ; i++) {
			// si pos[i] != 0, pos[i] = 1 car X^a = X
			if(pos[i] != 0) {
				res += P[i];
			}
		}

		return res;
	}

	for(i = 0 ; i < m ; i++) {
		pos_tmp = pos[i] % (p-1);
		if( (pos_tmp == 0) && (pos[i] != 0)) {
			pos_tmp = p-1;
		}
		res += pos_tmp*P[i];
	}

	return res;
}


// X_i^p = X_i
void nmod_multi_poly_mod_fermat(mp_ptr f_bar, const mp_ptr f, slong p, slong d, slong m, const slong* P) {

	if(p < d) {

		slong i, cpt, *pos, coord_p;
		
		_nmod_vec_zero(f_bar,n_pow(p,m));

		pos = (slong*)calloc(sizeof(slong),m);
		coord_p;
		cpt = 0;

		while(pos[m-1] != d) {

			coord_p = Coord_mod(P,p,pos,m);

			f_bar[coord_p] = n_addmod(f_bar[coord_p],f[cpt],p);

			pos[0]++;
			cpt++;
			i = 0;
			while( (pos[i] == d) && (i != (m-1)) ) {
				pos[i] = 0;
				// Le i est incrémenté avant pos++
				pos[++i]++;
			}
		}

		free(pos);

	}


}

