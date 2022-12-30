use std::{
    fmt::{self, Display},
    ops::{Add, AddAssign, Mul},
};

use rand::Rng;

use crate::N;
use crate::{TerPolyCoeffEncoding, HVC_MODULUS};

use super::{hvc::HVCPoly, hvc_ntt::HVCNTTPoly, lift};

#[derive(Debug, Clone, PartialEq, Copy)]
// A signed polynomial of degree N
pub struct SignedPoly {
    pub(crate) coeffs: [u32; N as usize],
}

impl Default for SignedPoly {
    fn default() -> Self {
        Self {
            coeffs: [0u32; N as usize],
        }
    }
}

impl Display for SignedPoly {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:#16x} {:#16x} {:#16x} {:#16x}",
            self.coeffs[0], self.coeffs[1], self.coeffs[3], self.coeffs[4]
        )
    }
}

impl From<&HVCPoly> for SignedPoly {
    fn from(input: &HVCPoly) -> Self {
        let mut res = Self::default();
        res.coeffs
            .iter_mut()
            .zip(input.coeffs)
            .for_each(|(e, f)| *e = f as i32);

        res
    }
}

impl From<&SignedPoly> for HVCPoly {
    fn from(input: &SignedPoly) -> Self {
        let mut res = Self::default();
        res.coeffs
            .iter_mut()
            .zip(input.coeffs)
            .for_each(|(e, f)| *e = lift(f, HVC_MODULUS) as u32);
        res
    }
}

impl From<&SignedPoly> for HVCNTTPoly {
    // convert poly into its ntt form. Requires that coefficients are between 0 and 12289
    fn from(poly: &SignedPoly) -> Self {
        (&HVCPoly::from(poly)).into()
    }
}

impl From<SignedPoly> for HVCPoly {
    fn from(input: SignedPoly) -> Self {
        (&input).into()
    }
}

impl From<HVCPoly> for SignedPoly {
    fn from(input: HVCPoly) -> Self {
        (&input).into()
    }
}

impl Add for SignedPoly {
    type Output = Self;

    // Coefficient wise additions without mod reduction.
    fn add(self, other: Self) -> Self {
        let mut res = Self::default();

        #[cfg(debug_assertions)]
        for (e, (f, g)) in res
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs))
        {
            *e = match i32::checked_add(*f, g) {
                Some(p) => p,
                None => {
                    panic!("overflowing additions")
                }
            };
        }

        #[cfg(not(debug_assertions))]
        for (e, (f, g)) in res
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *e = f + g;
        }
        res
    }
}

impl Mul for SignedPoly {
    type Output = Self;

    // Ring multiplication
    fn mul(self, other: Self) -> Self {
        (&(HVCPoly::from(&self) * HVCPoly::from(&other))).into()
    }
}

impl AddAssign for SignedPoly {
    // Coefficient wise additions without mod reduction.
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)]
        for (x, y) in self.coeffs.iter_mut().zip(other.coeffs) {
            *x = match i32::checked_add(*x, y) {
                Some(p) => p,
                None => {
                    panic!("overflowing additions")
                }
            };
        }
        #[cfg(not(debug_assertions))]
        for (x, y) in self.coeffs.iter_mut().zip(other.coeffs) {
            *x = *x + y
        }
    }
}

impl SignedPoly {
    // school book multiplication
    // slow. only used for correctness checking
    #[cfg(test)]
    pub(crate) fn schoolbook(a: &Self, b: &Self, q: i32) -> Self {
        let mut buf = [0i32; N as usize * 2];
        let mut c = [0; N as usize];
        for i in 0..N as usize {
            for j in 0..N as usize {
                buf[i + j] += (a.coeffs[i] as i64 * b.coeffs[j] as i64 % q as i64) as i32;
            }
        }
        for i in 0..N as usize {
            c[i] = (buf[i] - buf[i + N as usize]) % q;
        }
        Self { coeffs: c }
    }

    pub(crate) fn lifted<P: From<Self>>(&self, modulus: u32) -> P {
        let mut res = *self;
        res.normalize(modulus);
        res.into()
    }

    pub(crate) fn normalize(&mut self, modulus: u32) {
        for e in self.coeffs.iter_mut() {
            *e = lift(*e, modulus) as i32
        }
    }

    pub(crate) fn is_ternary(&self) -> bool {
        for &e in self.coeffs.iter() {
            if e != 0 && e != 1 && e != -1 {
                return false;
            }
        }
        true
    }

    // sample a random ternary polynomial with a fixed weight
    pub fn rand_ternary<R: Rng>(rng: &mut R, half_weight: usize) -> Self {
        let mut ct = 0;
        let mut coeffs = [0; N as usize];
        let mut rng_ct = 0;
        let mut tmp = rng.next_u32();

        while ct < half_weight {
            let index = (tmp & 0xFF) as usize;
            tmp >>= 9;
            rng_ct += 1;
            if rng_ct == 3 {
                tmp = rng.next_u32();
                rng_ct = 0;
            }
            if coeffs[index] == 0 {
                ct += 1;
                coeffs[index] = 1
            }
        }
        ct = 0;
        while ct < half_weight {
            let index = (tmp & 0xFF) as usize;
            tmp >>= 9;
            rng_ct += 1;
            if rng_ct == 3 {
                tmp = rng.next_u32();
                rng_ct = 0;
            }

            if coeffs[index] == 0 {
                ct += 1;
                coeffs[index] = -1
            }
        }
        Self { coeffs }
    }

    // sample a random binary polynomial
    pub fn rand_binary<R: Rng>(rng: &mut R) -> Self {
        let mut res = Self::default();
        for i in 0..16 {
            let mut tmp = rng.next_u32();
            for j in 0..32 {
                res.coeffs[i * 32 + j] = (tmp & 1) as i32;
                tmp >>= 1;
            }
        }

        res
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
