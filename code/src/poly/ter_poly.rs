#[derive(Debug, Clone, PartialEq)]
// ternary polynomials in coefficient encoding
pub struct TerPolyCoeffEncoding {
    pub(crate) pos: Vec<usize>,
    pub(crate) neg: Vec<usize>,
}

#[cfg(test)]
mod test {
    use super::TerPolyCoeffEncoding;
    use crate::{HOTSPoly, HVCPoly, Polynomial, ALPHA};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    #[test]
    fn test_ter_mul() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let half_weight = ALPHA / 2;

        for _ in 0..10 {
            let ter_poly = HVCPoly::rand_balanced_ternary(&mut rng, half_weight);
            let bin_poly = HVCPoly::rand_binary(&mut rng);
            let ter_poly_coeff_encoding: TerPolyCoeffEncoding = (&ter_poly).into();

            let prod_1 = HVCPoly::schoolbook(&bin_poly, &ter_poly);
            let prod_2 = HVCPoly::ter_mul_bin(&ter_poly_coeff_encoding, &bin_poly);
            let prod_3 = ter_poly * bin_poly;
            let prod_4 = HOTSPoly::from(&ter_poly_coeff_encoding) * HOTSPoly::from(&bin_poly);
            assert_eq!(prod_1, prod_3);
            assert_eq!(prod_1, prod_2);
            assert_eq!(HOTSPoly::from(&prod_1), prod_4)
        }
    }
}
