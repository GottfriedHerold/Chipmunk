use ark_std::{end_timer, start_timer};
use rand::Rng;

use crate::poly::Polynomial;
use crate::{HVCNTTPoly, HVCPoly, HVC_WIDTH};
#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct HVCHash {
    pub(crate) param_h: [HVCNTTPoly; HVC_WIDTH],
}

impl HVCHash {
    pub fn init<R: Rng>(rng: &mut R) -> Self {
        let timer = start_timer!(|| "initialize hvc hash");

        let mut res = Self::default();
        for e in res.param_h.iter_mut() {
            let tmp = HVCPoly::rand_poly(rng);
            *e = (&tmp).into();
        }

        end_timer!(timer);
        res
    }

    /// Hash function.
    /// Cost: 2*SMALL_MODULUS_BITS NTT and 1 INV_NTT.
    pub fn hash(&self, inputs: &[HVCPoly]) -> HVCPoly {
        // TODO: check the cost for fixed bases
        // may be faster than NTT
        assert_eq!(inputs.len(), HVC_WIDTH << 1);

        let mut res = HVCNTTPoly::default();
        let mut inputs = inputs.to_vec();
        inputs.iter_mut().for_each(|x| x.lift());

        #[cfg(feature = "parallel")]
        let prod_ntt: Vec<HVCNTTPoly> = self
            .param_h
            .iter()
            .zip(inputs.iter())
            .map(|(x, y)| (*x, *y))
            .collect::<Vec<(HVCNTTPoly, HVCPoly)>>()
            .into_par_iter()
            .map(|(x, y)| x * HVCNTTPoly::from(&y))
            .collect();

        #[cfg(not(feature = "parallel"))]
        let prod_ntt: Vec<HVCNTTPoly> = self
            .param_h
            .iter()
            .zip(inputs.iter())
            .map(|(x, y)| *x * HVCNTTPoly::from(y))
            .collect();

        for e in prod_ntt {
            res += e;
        }

        // convert the polynomial from NTT domain back to integers
        (&res).into()
    }

    pub(crate) fn decom_then_hash(&self, left: &HVCPoly, right: &HVCPoly) -> HVCPoly {
        self.hash_separate_inputs(&left.decompose(), &right.decompose())
    }

    pub(crate) fn hash_separate_inputs(&self, left: &[HVCPoly], right: &[HVCPoly]) -> HVCPoly {
        assert_eq!(left.len(), HVC_WIDTH);
        assert_eq!(right.len(), HVC_WIDTH);

        self.hash([left, right].concat().as_ref())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_hash() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

        let hasher = HVCHash::init(&mut rng);

        let inputs: Vec<HVCPoly> = (0..HVC_WIDTH << 1)
            .map(|_| HVCPoly::rand_poly(&mut rng).into())
            .collect();
        let _ = hasher.hash(&inputs);
    }

    #[test]
    fn test_homomorphism() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let hasher = HVCHash::init(&mut rng);

        for _ in 0..10 {
            // additive homomorphic
            {
                let poly1 = HVCPoly::rand_poly(&mut rng);
                let decomposed_poly1 = poly1.decompose();
                let poly2 = HVCPoly::rand_poly(&mut rng);
                let decomposed_poly2 = poly2.decompose();
                let r1 = HVCPoly::rand_ternary(&mut rng, 10);
                let r2 = HVCPoly::rand_ternary(&mut rng, 10);

                let decomposed_poly1: Vec<HVCPoly> =
                    decomposed_poly1.iter().map(|&x| r1 * x).collect();
                let decomposed_poly2: Vec<HVCPoly> =
                    decomposed_poly2.iter().map(|&x| r2 * x).collect();
                let decomposed: Vec<HVCPoly> = decomposed_poly1
                    .iter()
                    .zip(decomposed_poly2.iter())
                    .map(|(&x, &y)| x + y)
                    .collect();
                let poly_rec = HVCPoly::projection(&decomposed);
                let poly = poly1 * r1.into() + poly2 * r2.into();
                assert_eq!(poly, poly_rec);
            }
            // hash is homomorphic
            {
                let r1 = HVCPoly::rand_ternary(&mut rng, 10);
                let r2 = HVCPoly::rand_ternary(&mut rng, 10);

                let poly11 = HVCPoly::rand_poly(&mut rng);
                let poly12 = HVCPoly::rand_poly(&mut rng);
                let poly11_randomized: Vec<HVCPoly> =
                    poly11.decompose().iter().map(|&x| x * r1).collect();
                let poly12_randomized: Vec<HVCPoly> =
                    poly12.decompose().iter().map(|&x| x * r2).collect();

                let poly21 = HVCPoly::rand_poly(&mut rng);
                let poly22 = HVCPoly::rand_poly(&mut rng);
                let poly21_randomized: Vec<HVCPoly> =
                    poly21.decompose().iter().map(|&x| r2 * x).collect();
                let poly22_randomized: Vec<HVCPoly> =
                    poly22.decompose().iter().map(|&x| r2 * x).collect();

                let poly1 = hasher.decom_then_hash(&poly11, &poly12);
                let poly2 = hasher.decom_then_hash(&poly21, &poly22);
                let poly = poly1 * r1.into() + poly2 * r2.into();

                let polyx1_randomized: Vec<HVCPoly> = poly11_randomized
                    .iter()
                    .zip(poly21_randomized.iter())
                    .map(|(&x, &y)| x + y)
                    .collect();
                let polyx2_randomized: Vec<HVCPoly> = poly12_randomized
                    .iter()
                    .zip(poly22_randomized.iter())
                    .map(|(&x, &y)| x + y)
                    .collect();
                let poly_rec = hasher.hash_separate_inputs(&polyx1_randomized, &polyx2_randomized);

                assert_eq!(poly, poly_rec);
            }
        }
    }
}
