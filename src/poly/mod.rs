mod derive;
mod hots;
mod hvc;
mod ter_poly;

use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Mul},
};

pub use hots::HOTSNTTPoly;
pub use hots::HOTSPoly;
pub use hvc::HVCNTTPoly;
pub use hvc::HVCPoly;
use rand::Rng;
pub use ter_poly::TerPolyCoeffEncoding;

#[inline]
pub(crate) fn lift(a: i32, modulus: i32) -> i32 {
    (a % modulus + modulus) % modulus
}

#[inline]
pub(crate) fn normalize(a: i32, modulus: i32) -> i32 {
    let mut a = a % modulus;
    if a > modulus / 2 {
        a -= modulus
    }
    if a < -modulus / 2 {
        a += modulus
    }
    a
}

pub trait Polynomial:
    Sized + Add + Mul + Debug + Clone + Copy + Default + PartialEq + Display + AddAssign
{
    const DIM: usize;
    const MODULUS: i32;
    const MODULUS_OVER_2: i32;
    const SAMPLE_THRESHOLD: u32;

    type NTTPoly: From<Self> + Into<Self>;

    /// School book multiplication
    fn schoolbook(a: &Self, b: &Self) -> Self;

    /// sample a uniformly random polynomial with coefficients between 0 and q-1
    fn rand_poly<R: Rng>(rng: &mut R) -> Self;

    /// A 32 bytes digest of the polynomial
    fn digest(&self) -> [u8; 32];

    /// If the polynomial's coefficients are ternary
    fn is_ternary(&self) -> bool;

    /// Sample a random ternary polynomial with a fixed weight
    fn rand_balanced_ternary<R: Rng>(rng: &mut R, half_weight: usize) -> Self;

    /// Sample a random binary polynomial
    fn rand_binary<R: Rng>(rng: &mut R) -> Self;

    /// Sample a random ternary polynomial
    fn rand_ternary<R: Rng>(rng: &mut R, weight: usize) -> Self;

    /// Sample a random polynomial with coefficients between [-p, p]
    fn rand_mod_p<R: Rng>(rng: &mut R, p: u32) -> Self;

    /// Hash a blob into a message polynomial
    fn from_hash_message(msg: &[u8]) -> Self;

    /// Infinity norm of the polynomial
    fn infinity_norm(&self) -> u32;

    /// Normalize self into a polynomial within [-HVC_MODULUS_OVER_2, HVC_MODULUS_OVER_2)
    fn lift(&mut self);

    /// Normalize self into a polynomial within [-HVC_MODULUS_OVER_2, HVC_MODULUS_OVER_2)
    fn normalize(&mut self);
}

pub trait ChipmunkPolynomial: Polynomial {
    const WIDTH: usize;
    const ARITY: usize;

    /// Target polynomial that we are decomposing into.
    type TargetPoly;

    /// decompose a mod q polynomial into binary polynomials
    fn decompose(&self) -> Vec<Self::TargetPoly>;

    /// project a set of vectors to R_q
    fn projection(decomposed_polys: &[Self::TargetPoly]) -> Self;
}
