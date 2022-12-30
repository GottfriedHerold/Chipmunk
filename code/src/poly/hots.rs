use super::{lift, HOTSNTTPoly, normalize};
use crate::{
    HVCPoly, TerPolyCoeffEncoding, ALPHA_H, HOTS_MODULUS, HOTS_MODULUS_OVER_TWO,
    HOTS_SAMPLE_THRESHOLD, HOTS_WIDTH, N, PHI_ALPHA_H, TWO_ZETA_PLUS_ONE,
};
use crate::{HVC_MODULUS, HVC_MODULUS_OVER_TWO};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use sha2::Digest;
use sha2::Sha256;
use std::{
    fmt::{self, Display},
    ops::{Add, AddAssign, Mul},
};

#[derive(Debug, Clone, Copy)]
// HVC polynomials in canonical encoding
pub struct HOTSPoly {
    pub(crate) coeffs: [i32; N as usize],
}

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

impl PartialEq for HOTSPoly {
    fn eq(&self, rhs: &HOTSPoly) -> bool {
        self.coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(&x, &y)| lift(x, HOTS_MODULUS) == lift(y, HOTS_MODULUS))
            .fold(true, |acc, x| acc & x)
    }
}

impl Default for HOTSPoly {
    fn default() -> Self {
        Self {
            coeffs: [0i32; N as usize],
        }
    }
}

impl Display for HOTSPoly {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:#16x} {:#16x} {:#16x} {:#16x}",
            self.coeffs[0], self.coeffs[1], self.coeffs[0], self.coeffs[1]
        )
    }
}

impl Add for HOTSPoly {
    type Output = Self;

    // Coefficient wise additions without mod reduction.
    fn add(self, other: Self) -> Self {
        let mut res = self;
        res += other;
        res
    }
}

impl AddAssign for HOTSPoly {
    // Coefficient wise additions with mod reduction.
    fn add_assign(&mut self, other: Self) {
        self.coeffs
            .iter_mut()
            .zip(other.coeffs)
            .for_each(|(x, y)| *x = (*x + y) % HOTS_MODULUS as i32)
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
    // school book multiplication
    // slow. only used for correctness checking
    #[cfg(test)]
    pub(crate) fn schoolbook(a: &Self, b: &Self) -> Self {
        let mut buf = [0i32; N as usize * 2];
        let mut c = [0; N as usize];
        for i in 0..N as usize {
            for j in 0..N as usize {
                buf[i + j] +=
                    ((a.coeffs[i] as i64) * (b.coeffs[j] as i64) % HOTS_MODULUS as i64) as i32;
            }
        }
        for i in 0..N as usize {
            c[i] = lift(buf[i] - buf[i + N as usize], HOTS_MODULUS) as i32;
        }
        Self { coeffs: c }
    }

    /// sample a uniformly random polynomial with coefficients between 0 and q-1
    pub fn rand_poly<R: Rng>(rng: &mut R) -> Self {
        let mut res = Self::default();
        for e in res.coeffs.iter_mut() {
            let mut tmp = rng.next_u32();
            while tmp >= HOTS_SAMPLE_THRESHOLD {
                tmp = rng.next_u32();
            }
            *e = (tmp % HOTS_MODULUS as u32) as i32 - HOTS_MODULUS_OVER_TWO
        }
        res
    }

    /// decompose a mod q polynomial into binary polynomials
    pub fn decompose(&self) -> [HVCPoly; HOTS_WIDTH] {
        let mut res = [HVCPoly::default(); HOTS_WIDTH];
        let mut base_coeffs: Vec<_> = self.coeffs.iter().map(|&x| normalize(x, HOTS_MODULUS)).collect();
        for poly in res.iter_mut() {
            for (tar_coeff, cur_coeff) in (*poly).coeffs.iter_mut().zip(base_coeffs.iter_mut()) {
                *tar_coeff = *cur_coeff % TWO_ZETA_PLUS_ONE as i32;
                (*cur_coeff) /= TWO_ZETA_PLUS_ONE as i32;
            }
        }
        res
    }

    /// project a set of vectors to R_q
    pub fn projection(decomposed_polys: &[HVCPoly]) -> Self {
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

    /// Normalize self into a polynomial within [0, HOTS_MODULUS)
    pub fn lift(&mut self) {
        self.coeffs.iter_mut().for_each(|x| {
            *x = lift(*x, HOTS_MODULUS);
        });
    }

    /// A 256 digest of the polynomial
    pub(crate) fn _digest(&self) -> [u8; 32] {
        let mut inputs = Vec::new();
        for e in self.coeffs {
            inputs.push((e & 0xFF) as u8);
            inputs.push(((e >> 8) & 0xFF) as u8);
        }
        let mut hasher = Sha256::new();
        hasher.update(inputs);
        let result = hasher.finalize();
        result.into()
    }

    pub(crate) fn is_ternary(&self) -> bool {
        for &e in self.coeffs.iter() {
            if e != 1 && e != -1 {
                return false;
            }
        }
        true
    }

    // sample a random ternary polynomial with a fixed weight
    pub fn rand_balanced_ternary<R: Rng>(rng: &mut R, half_weight: usize) -> Self {
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

    // sample a random ternary polynomial
    pub fn rand_ternary<R: Rng>(rng: &mut R, weight: usize) -> Self {
        let mut res = Self::default();
        let mut ct = 0;
        // todo: improve sampling rates
        while ct < weight {
            let tmp = rng.next_u32();
            if res.coeffs[(tmp % N as u32) as usize] == 0 {
                ct += 1;
                if (tmp >> 9) & 1 == 1 {
                    res.coeffs[(tmp % N as u32) as usize] = 1;
                } else {
                    res.coeffs[(tmp % N as u32) as usize] = -1;
                }
            }
        }
        res
    }

    pub fn rand_mod_p<R: Rng>(rng: &mut R, p: u32) -> Self {
        // todo: improve sampling rates
        let mut res = Self::default();
        let modulus = (p * 2 + 1) as u32;
        let threshold = u32::MAX / modulus * modulus;
        for e in res.coeffs.iter_mut() {
            let mut tmp = rng.next_u32();
            while tmp > threshold {
                tmp = rng.next_u32();
            }

            *e = (tmp % modulus) as i32 - p as i32;
        }

        res
    }

    /// sample a random polynomial with coefficients between [-phi * alpha_H, phi * alpha_H]
    pub fn rand_mod_phi_alpha<R: Rng>(rng: &mut R) -> Self {
        Self::rand_mod_p(rng, PHI_ALPHA_H as u32)
    }

    /// hash a blob into a message polynomial
    pub(crate) fn from_hash_message(msg: &[u8]) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(msg);
        let seed = hasher.finalize().into();
        let mut rng = ChaCha20Rng::from_seed(seed);

        Self::rand_ternary(&mut rng, ALPHA_H)
    }

    pub fn infinity_norm(&self) -> u32 {
        let mut norm = 0;
        self.coeffs
            .iter()
            .for_each(|&x| norm = std::cmp::max(norm, std::cmp::max(x, -x)));
        norm as u32
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
    use rand::{RngCore, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use crate::{HOTSPoly, HOTS_MODULUS};

    #[test]
    fn test_normalization() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        for _ in 0..10 {
            let mut poly = HOTSPoly::rand_poly(&mut rng);
            let mut poly2 = poly.clone();

            poly.lift();
            poly2
                .coeffs
                .iter_mut()
                .for_each(|x| *x = *x + ((rng.next_u32() % 100) * HOTS_MODULUS as u32) as i32);
            poly2.lift();
            for (e, f) in poly.coeffs.iter().zip(poly2.coeffs.iter()) {
                assert_eq!(e, f)
            }
        }
    }

    #[test]
    fn test_decomposition() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

        // test decomposition is correct
        for _ in 0..10 {
            let poly = HOTSPoly::rand_poly(&mut rng);
            let decomposed = poly.decompose();
            let poly_rec = HOTSPoly::projection(&decomposed);
            assert_eq!(poly, poly_rec);
        }

        // test decomposition is homomorphic
        for _ in 0..10 {
            let poly1 = HOTSPoly::rand_poly(&mut rng);
            let decomposed_poly1 = poly1.decompose();
            let poly2 = HOTSPoly::rand_poly(&mut rng);
            let decomposed_poly2 = poly2.decompose();

            let decomposed: Vec<_> = decomposed_poly1
                .iter()
                .zip(decomposed_poly2.iter())
                .map(|(&x, &y)| x + y)
                .collect();
            let poly_rec = HOTSPoly::projection(&decomposed);
            let poly = poly1 + poly2;

            assert_eq!(poly, poly_rec);
        }
    }
}
