use std::io::{Read, Write};

use crate::{
    normalize, HVCPoly, Polynomial, ENCODING_NORM_BOUND, HVC_MODULUS, HVC_WIDTH, N,
    TWO_ZETA_PLUS_ONE, ZETA,
};

// TODO: alpha 1/2/3 has structs. we may use a better algorithm to compress them.
#[derive(Debug, Clone, Default)]
pub struct EncodedPoly {
    hint: HVCPoly,
    pub(crate) a_star: HVCPoly,
    pub(crate) alphas: Vec<HVCPoly>,
}

impl EncodedPoly {
    pub(crate) fn serialize<W: Write>(&self, mut _writer: W) {
        // // we use a simple byte repr for serialization
        // // there is a lot of zeros which can be removed for optimization
        // self.hint
        //     .coeffs
        //     .iter()
        //     .for_each(|coeff| writer.write_all(coeff.to_be_bytes().as_ref()).unwrap());
        // self.a_star
        //     .coeffs
        //     .iter()
        //     .for_each(|coeff| writer.write_all(coeff.to_be_bytes().as_ref()).unwrap());
        // self.alpha_1
        //     .coeffs
        //     .iter()
        //     .for_each(|coeff| writer.write_all(coeff.to_be_bytes().as_ref()).unwrap());

        // self.alpha_2
        //     .coeffs
        //     .iter()
        //     .for_each(|coeff| writer.write_all(coeff.to_be_bytes().as_ref()).unwrap());

        // self.alpha_3
        //     .coeffs
        //     .iter()
        //     .for_each(|coeff| writer.write_all(coeff.to_be_bytes().as_ref()).unwrap());
    }

    pub(crate) fn deserialize<R: Read>(mut _reader: R) -> Self {
        let mut _res = Self::default();
        _res
    }

    pub(crate) fn encode(polys: &[HVCPoly]) -> EncodedPoly {
        // proj_r is the hint over R_q
        let proj_r = HVCPoly::projection_r(polys);
        // proj_eta_kappa is over ZZ
        let proj_eta_kappa = HVCPoly::projection_zz(polys);
        let a_star = proj_eta_kappa
            .coeffs
            .iter()
            .zip(proj_r.coeffs.iter())
            .map(|(&a, &b)| {
                assert!((a - b) % HVC_MODULUS == 0);
                (a - b) / HVC_MODULUS
            })
            .collect::<Vec<_>>();
        let a_star = HVCPoly {
            coeffs: a_star.try_into().unwrap(),
        };

        let dec_eta_kappa = proj_eta_kappa.decompose_zz();

        // we do not really need to compute delta_v during encoding
        // let delta_v = polys
        //     .iter()
        //     .zip(dec_eta_kappa.iter())
        //     .map(|(&a, &b)| a - b)
        //     .collect::<Vec<_>>();

        let alphas = extract_alphas(polys, &dec_eta_kappa);

        EncodedPoly {
            hint: proj_r,
            a_star,
            alphas,
        }
    }

    pub(crate) fn decode(&self) -> [HVCPoly; HVC_WIDTH] {
        let h_double_prime = self
            .hint
            .coeffs
            .iter()
            .zip(self.a_star.coeffs.iter())
            .map(|(a, b)| a + HVC_MODULUS * b)
            .collect::<Vec<_>>();

        let h_double_prime = HVCPoly {
            coeffs: h_double_prime.try_into().unwrap(),
        };

        let delta_v = rebuild_delta_v(&self.alphas);

        let dec_h_double_prime = h_double_prime.decompose_zz();

        let mut res = [HVCPoly::default(); HVC_WIDTH];
        res.iter_mut()
            .zip(dec_h_double_prime.iter().zip(delta_v.iter()))
            .for_each(|(a, (b, c))| {
                *a = *b + *c;
            });
        res
    }

    // enforce norm bound for a_star and alphas
    pub(crate) fn is_norm_bounded(&self) -> bool {
        if self.a_star.infinity_norm() > ENCODING_NORM_BOUND {
            return false;
        }
        for alpha in self.alphas.iter() {
            if alpha.infinity_norm() > ENCODING_NORM_BOUND {
                return false;
            }
        }

        true
    }
}

// give the input vector v, and dec_{eta, kappa}(proj_{eta,kappa}(v))
// compute alphas s.t.
// delta_v = alpha1 b1 + alpha2 b2 + alpha3 b3
fn extract_alphas(polys: &[HVCPoly], dec_eta_kappa: &[HVCPoly]) -> Vec<HVCPoly> {
    let len = polys.len();
    let mut alphas = vec![HVCPoly::default(); len];

    for i in 0..N {
        // alpha_1 = -(v1 - w1 - a_star * q)/(2 * eta + 1)
        let tmp = polys[0].coeffs[i] - dec_eta_kappa[0].coeffs[i];
        assert!(tmp % TWO_ZETA_PLUS_ONE as i32 == 0);
        alphas[0].coeffs[i] = -tmp / TWO_ZETA_PLUS_ONE as i32;

        for j in 1..alphas.len()  {
            // alpha_2 = -(v2 - w2 - alpha_1)/(2 * eta + 1)
            let tmp = polys[j].coeffs[i] - dec_eta_kappa[j].coeffs[i] - alphas[j - 1].coeffs[i];
            assert!(tmp % TWO_ZETA_PLUS_ONE as i32 == 0);
            alphas[j].coeffs[i] = -tmp / TWO_ZETA_PLUS_ONE as i32;
        }

        // // alpha_3 = (v3 - w3 - alpha_2)/(2 * eta + 1)
        // let tmp =
        //     polys[len - 1].coeffs[i] - dec_eta_kappa[len - 1].coeffs[i] - alphas[len - 2].coeffs[i];
        // assert!(tmp % TWO_ZETA_PLUS_ONE as i32 == 0);
        // let r = tmp / TWO_ZETA_PLUS_ONE as i32;
        // alphas[len - 1].coeffs[i] = r;
    }
    alphas
}

// TODO: optimize this code
// given the alpha1/2/3
// compute delta_v s.t.
// delta_v = alpha1 b1 + alpha2 b2 + alpha3 b3
fn rebuild_delta_v(alphas: &[HVCPoly]) -> Vec<HVCPoly> {
    let len = alphas.len();
    let mut delta_vs = alphas.to_vec().clone();
    for delta_v in delta_vs.iter_mut() {
        for coeff in delta_v.coeffs.iter_mut() {
            *coeff *=-59
        }
    }

    for i in 0..len-1 {
        delta_vs[i+1] += alphas[i]
    }

    delta_vs
}

#[cfg(test)]
mod tests {
    use ark_std::{end_timer, start_timer};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use crate::Polynomial;

    use super::*;

    #[test]
    fn test_encoding() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let randomizer1 = HVCPoly::rand_balanced_ternary(&mut rng, 50);
        let randomizer2 = HVCPoly::rand_balanced_ternary(&mut rng, 50);
        let poly1 = HVCPoly::rand_poly(&mut rng);
        let poly2 = HVCPoly::rand_poly(&mut rng);
        let decomposed1 = poly1.decompose_r();
        let decomposed2 = poly2.decompose_r();
        let poly_rec = randomizer1 * poly1 + randomizer2 * poly2;

        let decomposed_sum = decomposed1
            .iter()
            .zip(decomposed2.iter())
            .map(|(&a, &b)| randomizer1 * a + randomizer2 * b)
            .collect::<Vec<_>>();

        let poly_rec1 = HVCPoly::projection_zz(&decomposed_sum);

        assert_eq!(poly_rec, poly_rec1);

        let encoded = EncodedPoly::encode(&decomposed_sum);
        let decoded = encoded.decode();
        assert_eq!(decomposed_sum, decoded);

        println!("encoded: {:?}", encoded);
    }

    #[test]
    fn bench_encoding() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

        let size = 1000;
        let mut decomposed_sums = vec![];

        for _ in 0..size {
            let randomizer1 = HVCPoly::rand_balanced_ternary(&mut rng, 10);
            let randomizer2 = HVCPoly::rand_balanced_ternary(&mut rng, 10);
            let poly1 = HVCPoly::rand_poly(&mut rng);
            let poly2 = HVCPoly::rand_poly(&mut rng);
            let decomposed1 = poly1.decompose_r();
            let decomposed2 = poly2.decompose_r();
            let poly_rec = randomizer1 * poly1 + randomizer2 * poly2;

            let decomposed_sum = decomposed1
                .iter()
                .zip(decomposed2.iter())
                .map(|(&a, &b)| randomizer1 * a + randomizer2 * b)
                .collect::<Vec<_>>();

            let poly_rec1 = HVCPoly::projection_zz(&decomposed_sum);

            assert_eq!(poly_rec, poly_rec1);
            decomposed_sums.push(decomposed_sum);
        }

        let mut encodeds = vec![];
        let encoding = start_timer!(|| format!("encode {} times", size));
        for i in 0..size {
            let encoded = EncodedPoly::encode(&decomposed_sums[i]);
            encodeds.push(encoded);
        }
        end_timer!(encoding);

        let encoding = start_timer!(|| format!("decode {} times", size));
        for i in 0..size {
            let decoded = encodeds[i].decode();
            assert_eq!(decomposed_sums[i], decoded);
        }
        end_timer!(encoding);
    }
}
