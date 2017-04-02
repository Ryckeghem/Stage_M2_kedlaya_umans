

mp_limb_t point_fixe(mp_limb_t mod, slong d, slong m, const mp_limb_t* vect_prime) {

	slong j, dm_pow = n_pow(d,m), dm = (d-1)*m+1;
	mp_limb_t point_fixe = mod;

	fmpz_t var, prod;
	fmpz_init(var);
	fmpz_init(prod);


	while(1) {

		fmpz_set_ui(var, point_fixe-1);
		fmpz_pow_ui(var, var, dm);
		fmpz_mul_ui(var, var, dm_pow);

		// Calcul de k
		fmpz_set_ui(prod, 1);
		j = 0;

		while(fmpz_cmp(prod, var) < 0) {
			fmpz_mul_ui(prod, prod, vect_prime[j]);
			j++;
		}

		if(point_fixe <= vect_prime[j-1]) {
			fmpz_clear(var);
			fmpz_clear(prod);
			return point_fixe;
		}

		point_fixe = vect_prime[j-1];

	}

}
