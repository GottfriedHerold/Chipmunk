use crate::normalize;
use crate::{impl_ntt_poly, impl_poly, HOTS_ONE_OVER_N};
use crate::{HOTSPoly, HOTS_MODULUS, N};

use super::ntt_param::{INV_NTT_TABLE, NTT_TABLE};

#[derive(Debug, Clone, Copy)]
// HVC polynomials in NTT encoding
pub struct HOTSNTTPoly {
    pub(crate) coeffs: [i32; N as usize],
}

impl_poly!(HOTSNTTPoly, HOTS_MODULUS, N);
impl_ntt_poly!(HOTSPoly, HOTSNTTPoly, HOTS_MODULUS, N, HOTS_ONE_OVER_N);

#[cfg(test)]
mod test {
    use super::*;
    use super::{inv_ntt, ntt};
    use crate::impl_ntt_poly_tests;
    use crate::Polynomial;
    use crate::{poly::lift, HOTSPoly, HOTS_MODULUS};
    use ark_std::{end_timer, start_timer};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    impl_ntt_poly_tests!(HOTSPoly, HOTSNTTPoly, HOTS_MODULUS);
}
