

void nmod_multi_poly_clear(nmod_multi_poly_t poly) {
    if (poly->coeffs)
        flint_free(poly->coeffs);
}
