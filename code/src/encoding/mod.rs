use std::io::{Read, Write};

use crate::{normalize, HVCPoly, Polynomial, HVC_MODULUS, HVC_WIDTH, TWO_ZETA_PLUS_ONE, ZETA};

// TODO: alpha 1/2/3 has structs. we may use a better algorithm to compress them.
#[derive(Debug, Clone, Default)]
pub struct EncodedPoly {
    // todo: delta_v, can be reconstructed from the rest
    delta_v: Vec<HVCPoly>,
    hint: HVCPoly,
    a_star: Vec<i32>,
    alpha_1: Vec<i32>,
    alpha_2: Vec<i32>,
    alpha_3: Vec<i32>,
}

impl EncodedPoly {
    pub(crate) fn serialize<W: Write>(&self, mut _writer: W) {
        // not implemented yet
    }

    pub(crate) fn deserialize<R: Read>(mut _reader: R) -> Self {
        // not implemented yet
        Self::default()
    }

    pub(crate) fn encode(polys: &[HVCPoly]) -> EncodedPoly {
        // proj_r is the hint over R_q
        let proj_r = HVCPoly::projection_r(polys);
        // let mut hint = proj_r;
        // hint.normalize();
        // proj_eta_kappa is over ZZ
        let proj_eta_kappa = HVCPoly::projection_zz(polys);
        let a_star = proj_eta_kappa
            .coeffs
            .iter()
            .zip(proj_r.coeffs.iter())
            .map(|(&a, &b)| {
                // comment this line out for optimization
                assert!((a - b) % HVC_MODULUS == 0);
                normalize((a - b) / HVC_MODULUS, TWO_ZETA_PLUS_ONE as i32)
            })
            .collect::<Vec<_>>();

        let dec_eta_kappa = proj_eta_kappa.decompose_zz();

        let delta_v = polys
            .iter()
            .zip(dec_eta_kappa.iter())
            .map(|(&a, &b)| a - b)
            .collect::<Vec<_>>();

        // println!("proj_eta_kappa: {:?}", proj_eta_kappa);
        // println!("proj_r: {:?}", proj_r);
        // println!("hint': {:?}", hint);
        // println!("a star: {:?}", a_star);
        // println!("dec_eta_kappa: {:?}", dec_eta_kappa);
        // println!("delta_v: {:?}", delta_v);

        // println!("v1: {:?}", polys[0]);
        // println!("w1: {:?}", dec_eta_kappa[0]);

        // alpha_1 = -(v1 - w1 - a_star * q)/(2 * eta + 1)
        let mut alpha_1 = vec![];
        for i in 0..512 {
            // let tmp = normalize(polys[0].coeffs[i] - dec_eta_kappa[0].coeffs[i], HVC_MODULUS);
            let tmp = polys[0].coeffs[i] - dec_eta_kappa[0].coeffs[i];
            // println!(
            //     "{} {} {} {} {} {}",
            //     i,
            //     polys[0].coeffs[i],
            //     a_star[i],
            //     dec_eta_kappa[0].coeffs[i],
            //     tmp,
            //     tmp % 59
            // );

            assert!(tmp % TWO_ZETA_PLUS_ONE as i32 == 0);
            alpha_1.push(-tmp / TWO_ZETA_PLUS_ONE as i32);
        }

        // println!("alpha1: {:?}", alpha_1);

        // alpha_2 = -(v2 - w2 - alpha_1)/(2 * eta + 1)
        let mut alpha_2 = vec![];
        for i in 0..512 {
            // let tmp = normalize(
            //     polys[1].coeffs[i] - dec_eta_kappa[1].coeffs[i] - alpha_1[i],
            //     HVC_MODULUS,
            // );
            let tmp = polys[1].coeffs[i] - dec_eta_kappa[1].coeffs[i] - alpha_1[i];
            // println!(
            //     "{} {} {} {} {} {}",
            //     i, polys[1].coeffs[i], a_star[i], dec_eta_kappa[1].coeffs[i], alpha_1[i], tmp
            // );
            assert!(tmp % TWO_ZETA_PLUS_ONE as i32 == 0);
            alpha_2.push(-tmp / TWO_ZETA_PLUS_ONE as i32);
        }

        // println!("alpha2: {:?}", alpha_2);

        // alpha_3 = (v3 - w3 - alpha_2)/(2 * eta + 1)
        let mut alpha_3 = vec![];
        for i in 0..512 {
            // let tmp = normalize(
            //     polys[2].coeffs[i] - dec_eta_kappa[2].coeffs[i] - alpha_2[i],
            //     HVC_MODULUS,
            // );
            let tmp = polys[2].coeffs[i] - dec_eta_kappa[2].coeffs[i] - alpha_2[i];
            let r = tmp / TWO_ZETA_PLUS_ONE as i32;

            alpha_3.push(r);
        }
        // println!("alpha3: {:?}", alpha_3);

        EncodedPoly {
            hint: proj_r,
            delta_v,
            a_star,
            alpha_1,
            alpha_2,
            alpha_3,
        }
    }

    pub(crate) fn decoding(&self) -> [HVCPoly; HVC_WIDTH] {
        let h_double_prime = self
            .hint
            .coeffs
            .iter()
            .zip(self.a_star.iter())
            .map(|(a, b)| a + HVC_MODULUS * b)
            .collect::<Vec<_>>();
        let dec_h_double_prime = HVCPoly {
            coeffs: h_double_prime.try_into().unwrap(),
        }
        .decompose_zz();
        let mut res = [HVCPoly::default(); HVC_WIDTH];

        res.iter_mut()
            .zip(dec_h_double_prime.iter())
            .for_each(|(a, b)| *a += *b);
        res
    }
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
        let decomposed1 = poly1.decompose_r();
        let decomposed2 = poly2.decompose_r();
        let decomposed_sum = decomposed1
            .iter()
            .zip(decomposed2.iter())
            .map(|(&a, &b)| randomizer1 * a + randomizer2 * b)
            // .map(|(&a, &b)| a+ b)
            .collect::<Vec<_>>();

        let poly_rec = HVCPoly::projection_r(&decomposed_sum);
        let encoded = EncodedPoly::encode(&decomposed_sum);
        println!("encoded: {:?}", encoded);
        let decoded = encoded.decoding();
        println!("works: {}", decomposed_sum == decoded.to_vec());
        for i in 0..512 {
            println!(
                "{} {} {} {}",
                i,
                decomposed_sum[0].coeffs[i],
                decoded[0].coeffs[i],
                (decomposed_sum[0].coeffs[i] - decoded[0].coeffs[i]) % 59
            );
        }

        println!("works: {}", poly_rec == HVCPoly::projection_r(&decoded));
    }
}
