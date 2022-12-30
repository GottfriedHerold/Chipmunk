mod hots;
mod hots_ntt;
mod hots_ntt_core;
mod hvc;
mod hvc_ntt;
mod ter_poly;

mod hvc_ntt_core;

pub use hots::HOTSPoly;
pub use hots_ntt::HOTSNTTPoly;
pub use hvc::HVCPoly;
pub use hvc_ntt::HVCNTTPoly;
pub use ter_poly::TerPolyCoeffEncoding;

pub(crate) use hots_ntt_core::*;
pub(crate) use hvc_ntt_core::*;

#[inline]
fn lift(a: i32, modulus: i32) -> i32 {
    (a % modulus + modulus) % modulus
}

#[inline]
fn normalize(a: i32, modulus: i32) -> i32 {
    let mut a = a % modulus;
    if a > modulus / 2 {
        a -= modulus
    }
    if a < -modulus / 2 {
        a += modulus
    }
    a
}
