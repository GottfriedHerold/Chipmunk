use super::lift;
use super::normalize;
use crate::impl_signed_poly;
use crate::impl_signed_poly_functions;
use crate::Polynomial;
use crate::{
    TerPolyCoeffEncoding, HVC_MODULUS, HVC_MODULUS_OVER_TWO, HVC_SAMPLE_THRESHOLD, HVC_WIDTH, N,
    TWO_ZETA_PLUS_ONE,
};
use rand::Rng;
use sha2::Digest;
use sha2::Sha256;
use std::ops::Mul;

mod hvc_ntt;
mod hvc_ntt_param;

pub use hvc_ntt::HVCNTTPoly;

#[derive(Debug, Clone, Copy)]
// HVC polynomials in canonical encoding
pub struct HVCPoly {
    pub(crate) coeffs: [i32; N as usize],
}

impl_signed_poly!(HVCPoly, HVC_MODULUS, N);
impl_signed_poly_functions!(
    HVCPoly,
    HVCNTTPoly,
    HVC_MODULUS,
    N,
    HVC_MODULUS_OVER_TWO,
    HVC_SAMPLE_THRESHOLD
);

impl Mul for HVCPoly {
    type Output = Self;

    // Ring multiplication
    fn mul(self, other: Self) -> Self {
        (&(HVCNTTPoly::from(&self) * HVCNTTPoly::from(&other))).into()
    }
}

impl HVCPoly {
    /// decompose a polynomial into binary polynomials
    pub fn decompose_zz(&self) -> [HVCPoly; HVC_WIDTH] {
        let mut res = [HVCPoly::default(); HVC_WIDTH];
        let mut base_coeffs: Vec<_> = self
            .coeffs
            .iter()
            .map(|&x| normalize(x, HVC_MODULUS))
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
        let mut res = decomposed_polys[HVC_WIDTH - 1];
        res.coeffs
            .iter_mut()
            .for_each(|x| *x = normalize(*x, HVC_MODULUS));
        for decomposed_poly in decomposed_polys.iter().rev().skip(1) {
            for (res, &base) in res.coeffs.iter_mut().zip(decomposed_poly.coeffs.iter()) {
                *res *= TWO_ZETA_PLUS_ONE as i32;
                *res += normalize(base, HVC_MODULUS);
            }
        }
        res
    }

    /// decompose a mod q polynomial into binary polynomials
    pub fn decompose_r(&self) -> [HVCPoly; HVC_WIDTH] {
        let mut res = [HVCPoly::default(); HVC_WIDTH];
        let mut base_coeffs: Vec<_> = self
            .coeffs
            .iter()
            .map(|&x| normalize(x, HVC_MODULUS))
            .collect();
        for poly in res.iter_mut() {
            for (tar_coeff, cur_coeff) in (*poly).coeffs.iter_mut().zip(base_coeffs.iter_mut()) {
                *tar_coeff = *cur_coeff % TWO_ZETA_PLUS_ONE as i32;
                (*cur_coeff) /= TWO_ZETA_PLUS_ONE as i32;
            }
        }
        res
    }

    /// project a set of vectors to R
    pub fn projection_r(decomposed_polys: &[HVCPoly]) -> Self {
        let mut res = decomposed_polys[HVC_WIDTH - 1];
        for decomposed_poly in decomposed_polys.iter().rev().skip(1) {
            for (res, &base) in res.coeffs.iter_mut().zip(decomposed_poly.coeffs.iter()) {
                *res *= TWO_ZETA_PLUS_ONE as i32;
                *res += base;
            }
        }
        res
    }

    /// Normalize self into a polynomial within [-HVC_MODULUS_OVER_2, HVC_MODULUS_OVER_2)
    pub fn lift(&mut self) {
        self.coeffs.iter_mut().for_each(|x| {
            *x = lift(*x, HVC_MODULUS);
        });
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
    use crate::poly::Polynomial;
    use rand::{RngCore, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use crate::{impl_poly_tests, HVCPoly, HVC_MODULUS};

    impl_poly_tests!(HVCPoly, HVC_MODULUS);
}
