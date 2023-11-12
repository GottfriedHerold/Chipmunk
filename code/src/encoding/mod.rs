use crate::{HVCPoly, HVC_MODULUS, HVC_WIDTH};

pub(crate) fn encode(polys: &[HVCPoly]) {
    // proj_r is the hint over R_q
    let proj_r = HVCPoly::projection(polys);
    // proj_eta_kappa is over R
    let proj_eta_kappa = HVCPoly::projection_r(polys);
    let a_star = proj_r
        .coeffs
        .iter()
        .zip(proj_eta_kappa.coeffs.iter())
        .map(|(&a, &b)| {
            // comment this line out for optimization
            assert!((a - b) % HVC_MODULUS == 0);
            (a - b) / HVC_MODULUS
        })
        .collect::<Vec<_>>();

    println!("proj_eta_kappa: {:?}", proj_eta_kappa);
    println!("proj_r: {:?}", proj_r);
    println!("diff: {:?}", a_star);
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use crate::{Polynomial, TerPolyCoeffEncoding};

    use super::*;

    #[test]
    fn test_encoding() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let randomizer1 = HVCPoly::rand_balanced_ternary(&mut rng, 10);
        let randomizer2 = HVCPoly::rand_balanced_ternary(&mut rng, 10);
        let poly1 = HVCPoly::rand_poly(&mut rng);
        let poly2 = HVCPoly::rand_poly(&mut rng);
        let decomposed1 = poly1.decompose();
        let decomposed2 = poly2.decompose();
        let decomposed_sum = decomposed1
            .iter()
            .zip(decomposed2.iter())
            .map(|(&a, &b)| randomizer1 * a + randomizer2 * b)
            // .map(|(&a, &b)| a+ b)
            .collect::<Vec<_>>();

        let poly_rec = HVCPoly::projection(&decomposed_sum);
        encode(&decomposed_sum);
    }
}
