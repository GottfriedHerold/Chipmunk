use super::{lift, normalize};
use crate::Polynomial;
use crate::{
    impl_signed_poly, impl_signed_poly_functions, HVCPoly, TerPolyCoeffEncoding, HOTS_MODULUS,
    HOTS_MODULUS_OVER_TWO, HOTS_SAMPLE_THRESHOLD, HOTS_WIDTH, N, PHI_ALPHA_H, TWO_ZETA_PLUS_ONE,
};
use crate::{HVC_MODULUS, HVC_MODULUS_OVER_TWO};
use rand::Rng;
use sha2::Digest;
use sha2::Sha256;
use std::ops::Mul;

mod hots_ntt;
mod ntt_param;

pub use hots_ntt::HOTSNTTPoly;

#[derive(Debug, Clone, Copy)]
// HVC polynomials in canonical encoding
pub struct HOTSPoly {
    pub(crate) coeffs: [i32; N as usize],
}

impl_signed_poly!(HOTSPoly, HOTS_MODULUS, N);
impl_signed_poly_functions!(
    HOTSPoly,
    HOTSNTTPoly,
    HOTS_MODULUS,
    N,
    HOTS_MODULUS_OVER_TWO,
    HOTS_SAMPLE_THRESHOLD
);

impl From<&HVCPoly> for HOTSPoly {
    fn from(p: &HVCPoly) -> Self {
        let mut coeffs = p.coeffs.clone();
        for x in coeffs.iter_mut() {
            *x = *x % HVC_MODULUS as i32;
            if *x > HVC_MODULUS_OVER_TWO {
                *x -= HVC_MODULUS as i32
            }
            if *x < -HVC_MODULUS_OVER_TWO {
                *x += HVC_MODULUS as i32
            }
        }

        Self { coeffs }
    }
}

impl Mul for HOTSPoly {
    type Output = Self;

    // Ring multiplication
    fn mul(self, other: Self) -> Self {
        (&(HOTSNTTPoly::from(&self) * HOTSNTTPoly::from(&other))).into()
    }
}

impl HOTSPoly {
    /// decompose a mod q polynomial into binary polynomials
    pub fn decompose_zz(&self) -> [HVCPoly; HOTS_WIDTH] {
        let mut res = [HVCPoly::default(); HOTS_WIDTH];
        let mut base_coeffs: Vec<_> = self
            .coeffs
            .iter()
            .map(|&x| normalize(x, HOTS_MODULUS))
            .collect();
        for poly in res.iter_mut() {
            for (tar_coeff, cur_coeff) in (*poly).coeffs.iter_mut().zip(base_coeffs.iter_mut()) {
                *tar_coeff = *cur_coeff % TWO_ZETA_PLUS_ONE as i32;
                (*cur_coeff) /= TWO_ZETA_PLUS_ONE as i32;
            }
        }
        res
    }

    /// project a set of vectors to R_q
    pub fn projection_zz(decomposed_polys: &[HVCPoly]) -> Self {
        let mut res = decomposed_polys[HOTS_WIDTH - 1];
        res.coeffs
            .iter_mut()
            .for_each(|x| *x = normalize(*x, HVC_MODULUS));
        for decomposed_poly in decomposed_polys.iter().rev().skip(1) {
            for (res, &base) in res.coeffs.iter_mut().zip(decomposed_poly.coeffs.iter()) {
                *res *= TWO_ZETA_PLUS_ONE as i32;
                *res += normalize(base, HVC_MODULUS);
            }
        }
        Self { coeffs: res.coeffs }
    }

    /// decompose a mod q polynomial into binary polynomials
    pub fn decompose_r(&self) -> [HVCPoly; HOTS_WIDTH] {
        let mut res = [HVCPoly::default(); HOTS_WIDTH];
        let mut base_coeffs: Vec<_> = self
            .coeffs
            .iter()
            .map(|&x| normalize(x, HOTS_MODULUS))
            .collect();
        for poly in res.iter_mut() {
            for (tar_coeff, cur_coeff) in (*poly).coeffs.iter_mut().zip(base_coeffs.iter_mut()) {
                *tar_coeff = *cur_coeff % TWO_ZETA_PLUS_ONE as i32;
                (*cur_coeff) /= TWO_ZETA_PLUS_ONE as i32;
            }
        }
        res
    }

    /// project a set of vectors to R_q
    pub fn projection_r(decomposed_polys: &[HVCPoly]) -> Self {
        let mut res = decomposed_polys[HOTS_WIDTH - 1];
        res.coeffs
            .iter_mut()
            .for_each(|x| *x = normalize(*x, HVC_MODULUS));
        for decomposed_poly in decomposed_polys.iter().rev().skip(1) {
            for (res, &base) in res.coeffs.iter_mut().zip(decomposed_poly.coeffs.iter()) {
                *res *= TWO_ZETA_PLUS_ONE as i32;
                *res += normalize(base, HVC_MODULUS);
            }
        }
        Self { coeffs: res.coeffs }
    }

    /// sample a random polynomial with coefficients between [-phi * alpha_H, phi * alpha_H]
    pub fn rand_mod_phi_alpha<R: Rng>(rng: &mut R) -> Self {
        Self::rand_mod_p(rng, PHI_ALPHA_H as u32)
    }

    // multiply a ternary with a binary poly
    pub fn ter_mul_bin(ter: &TerPolyCoeffEncoding, bin: &Self) -> Self {
        Self::from(ter) * *bin

        // TODO: use the following AVX code
        // #[cfg(debug_assertions)]
        // assert!(bin.is_binary());

        // let mut res = Self::default();
        // let mut tmp = [0i8; N];
        // let mut buf = [0u8; 2 * N];
        // let bin: Vec<i8> = bin.coeffs.iter().map(|&x| x as i8).collect();
        // let ter: Vec<u8> = ter.indices.iter().map(|&x| x as u8).collect();

        // unsafe {
        //     ternary_mul(
        //         tmp.as_mut_ptr(),
        //         buf.as_mut_ptr(),
        //         bin.as_ptr(),
        //         ter.as_ptr(),
        //     );
        // }
        // for (e, f) in res.coeffs.iter_mut().zip(tmp.iter()) {
        //     *e = *f as i32
        // }
        // res
    }
}

#[cfg(test)]
mod test {
    use crate::Polynomial;
    use rand::{RngCore, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use crate::{impl_poly_tests, HOTSPoly, HOTS_MODULUS};

    impl_poly_tests!(HOTSPoly, HOTS_MODULUS);
}
