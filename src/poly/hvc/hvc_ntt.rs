use super::hvc_ntt_param::{INV_NTT_TABLE, NTT_TABLE};
use crate::normalize;
use crate::{impl_ntt_poly, impl_poly, HVCPoly, HVC_ONE_OVER_N};
use crate::{HVC_MODULUS, N};

#[derive(Debug, Clone, Copy)]
// HVC polynomials in NTT encoding
pub struct HVCNTTPoly {
    pub(crate) coeffs: [i32; N as usize],
}

impl_poly!(HVCNTTPoly, HVC_MODULUS, N);
impl_ntt_poly!(HVCPoly, HVCNTTPoly, HVC_MODULUS, N, HVC_ONE_OVER_N);

#[cfg(test)]
mod test {
    use super::*;
    use super::{inv_ntt, ntt};
    use crate::impl_ntt_poly_tests;
    use crate::poly::Polynomial;
    use crate::{poly::lift, HVCPoly, HVC_MODULUS};
    use ark_std::{end_timer, start_timer};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    impl_ntt_poly_tests!(HVCPoly, HVCNTTPoly, HVC_MODULUS);
}
