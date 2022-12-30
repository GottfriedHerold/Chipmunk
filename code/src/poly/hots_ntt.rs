use std::ops::{Add, AddAssign, Mul};

use crate::{hots_inv_ntt, hots_ntt};

use crate::{HOTSPoly, HOTS_MODULUS, N};

#[derive(Debug, Clone, PartialEq, Copy)]
// HVC polynomials in NTT encoding
pub struct HOTSNTTPoly {
    pub(crate) coeffs: [i32; N as usize],
}

impl Default for HOTSNTTPoly {
    fn default() -> Self {
        Self {
            coeffs: [0i32; N as usize],
        }
    }
}

impl From<&HOTSPoly> for HOTSNTTPoly {
    // convert poly into its ntt form. Requires that coefficients are between 0 and 12289
    fn from(poly: &HOTSPoly) -> Self {
        let mut coeffs = poly.coeffs;
        hots_ntt(&mut coeffs);
        Self { coeffs }
    }
}

impl From<&HOTSNTTPoly> for HOTSPoly {
    fn from(poly: &HOTSNTTPoly) -> Self {
        let mut coeffs = poly.coeffs;
        hots_inv_ntt(&mut coeffs);
        HOTSPoly { coeffs }
    }
}

impl Add for HOTSNTTPoly {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut res = Self::default();
        for (e, (f, g)) in res
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *e = (f + g) % HOTS_MODULUS as i32
        }

        res
    }
}

impl AddAssign for HOTSNTTPoly {
    fn add_assign(&mut self, other: HOTSNTTPoly) {
        for (x, y) in self.coeffs.iter_mut().zip(other.coeffs) {
            *x = (*x + y) % HOTS_MODULUS as i32
        }
    }
}

impl Mul for HOTSNTTPoly {
    type Output = Self;

    // Coefficient-wise multiplication over the NTT domain.
    fn mul(self, other: Self) -> Self {
        let mut res = Self::default();
        for (e, (f, g)) in res
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *e = (((*f as i64) * (*g as i64)) % HOTS_MODULUS as i64) as i32
        }

        res
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_conversion() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        for _ in 0..10 {
            let poly = HOTSPoly::rand_poly(&mut rng);
            let poly_ntt: HOTSNTTPoly = (&poly).into();
            let poly_rec: HOTSPoly = (&poly_ntt).into();

            assert_eq!(poly, poly_rec)
        }
    }

    #[test]
    fn test_arithmetic() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        for _ in 0..10 {
            let a = HOTSPoly::rand_poly(&mut rng);
            let a_ntt: HOTSNTTPoly = (&a).into();
            let b = HOTSPoly::rand_poly(&mut rng);
            let b_ntt: HOTSNTTPoly = (&b).into();

            {
                // test correctness of ntt multiplications
                let c_ntt = a_ntt * b_ntt;
                let c: HOTSPoly = (&c_ntt).into();
                let c_rec = HOTSPoly::schoolbook(&a, &b);

                assert_eq!(c, c_rec);
            }
            {
                // test correctness of ntt additions
                let d_ntt = a_ntt + b_ntt;
                let d: HOTSPoly = (&d_ntt).into();
                let d_rec = a + b;

                assert_eq!(d, d_rec)
            }
        }
    }
}
