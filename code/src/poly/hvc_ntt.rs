use super::hvc::HVCPoly;
use crate::{hvc_inv_ntt, hvc_ntt};
use crate::{HVC_MODULUS, N};
use std::ops::{Add, AddAssign, Mul};

#[derive(Debug, Clone, PartialEq, Copy)]
// HVC polynomials in NTT encoding
pub struct HVCNTTPoly {
    pub(crate) coeffs: [i32; N as usize],
}

impl Default for HVCNTTPoly {
    fn default() -> Self {
        Self {
            coeffs: [0i32; N as usize],
        }
    }
}

impl From<&HVCPoly> for HVCNTTPoly {
    // convert poly into its ntt form. Requires that coefficients are between 0 and 12289
    fn from(poly: &HVCPoly) -> Self {
        let mut coeffs = poly.coeffs;
        hvc_ntt(&mut coeffs);
        Self { coeffs }
    }
}

impl From<&HVCNTTPoly> for HVCPoly {
    fn from(poly: &HVCNTTPoly) -> Self {
        let mut coeffs = poly.coeffs;
        hvc_inv_ntt(&mut coeffs);
        HVCPoly { coeffs }
    }
}

impl Add for HVCNTTPoly {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut res = Self::default();
        for (e, (f, g)) in res
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *e = (f + g) % HVC_MODULUS as i32
        }

        res
    }
}

impl AddAssign for HVCNTTPoly {
    fn add_assign(&mut self, other: HVCNTTPoly) {
        for (x, y) in self.coeffs.iter_mut().zip(other.coeffs) {
            *x = (*x + y) % HVC_MODULUS as i32
        }
    }
}

impl Mul for HVCNTTPoly {
    type Output = Self;

    // Coefficient-wise multiplication over the NTT domain.
    fn mul(self, other: Self) -> Self {
        let mut res = Self::default();
        for (e, (f, g)) in res
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *e = (((*f as i64) * (*g as i64)) % HVC_MODULUS as i64) as i32
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
            let poly = HVCPoly::rand_poly(&mut rng);
            let poly_ntt: HVCNTTPoly = (&poly).into();
            let poly_rec: HVCPoly = (&poly_ntt).into();

            assert_eq!(poly, poly_rec)
        }
    }

    #[test]
    fn test_arithmetic() {
        {
            let a = HVCPoly {
                coeffs: [1; N as usize],
            };
            let a_ntt: HVCNTTPoly = (&a).into();
            let b = HVCPoly {
                coeffs: [1; N as usize],
            };
            let b_ntt: HVCNTTPoly = (&b).into();

            // test correctness of ntt multiplications
            let c_ntt = a_ntt * b_ntt;
            let c: HVCPoly = (&c_ntt).into();
            let c_rec = HVCPoly::schoolbook(&a, &b);

            assert_eq!(c, c_rec);
        }

        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        for _ in 0..10 {
            let a = HVCPoly::rand_poly(&mut rng);
            let a_ntt: HVCNTTPoly = (&a).into();
            let b = HVCPoly::rand_poly(&mut rng);
            let b_ntt: HVCNTTPoly = (&b).into();

            {
                // test correctness of ntt multiplications
                let c_ntt = a_ntt * b_ntt;
                let c: HVCPoly = (&c_ntt).into();
                let c_rec = HVCPoly::schoolbook(&a, &b);

                assert_eq!(c, c_rec);
            }
            {
                // test correctness of ntt additions
                let d_ntt = a_ntt + b_ntt;
                let d: HVCPoly = (&d_ntt).into();
                let d_rec = a + b;

                assert_eq!(d, d_rec)
            }
        }
    }
}
