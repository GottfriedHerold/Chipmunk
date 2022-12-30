#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::ops::AddAssign;

use crate::{poly::HOTSPoly, randomizer::Randomizers, TerPolyCoeffEncoding, GAMMA};

// HOTS signature
#[derive(Debug, Clone, Copy)]
pub struct HotsSig {
    pub(crate) sigma: [HOTSPoly; GAMMA],
    pub(crate) is_randomized: bool,
}

impl HotsSig {
    /// Randomize an Hots Signature
    pub fn randomize_with(&mut self, ternary: &TerPolyCoeffEncoding) {
        if self.is_randomized {
            panic!("already randomized")
        }

        #[cfg(feature = "parallel")]
        self.sigma.par_iter_mut().for_each(|x| {
            *x = *x * HOTSPoly::from(ternary);
        });

        #[cfg(not(feature = "parallel"))]
        self.sigma.iter_mut().for_each(|x| {
            *x = *x * *ternary;
        });

        self.is_randomized = true;
    }

    /// aggregated randomized signatures
    pub(crate) fn aggregate_randomized_signatures(sigs: &[Self]) -> Self {
        let mut res = sigs[0];
        for &e in sigs.iter().skip(1) {
            res += e;
        }
        res
    }

    ///
    pub(crate) fn aggregate_with_randomizers(sigs: &[Self], randomizers: &Randomizers) -> Self {
        #[cfg(feature = "parallel")]
        {
            let mut sig_and_randomizers: Vec<(Self, TerPolyCoeffEncoding)> = sigs
                .iter()
                .zip(randomizers.poly.iter())
                .map(|(&s, r)| (s, r.clone()))
                .collect();

            sig_and_randomizers
                .iter_mut()
                .for_each(|(s, randomizer)| s.randomize_with(randomizer));
            let sig_randomized: Vec<HotsSig> =
                sig_and_randomizers.iter().map(|(s, _r)| *s).collect();
            Self::aggregate_randomized_signatures(&sig_randomized)
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut sig_randomized = sigs.to_vec();

            sig_randomized
                .iter_mut()
                .zip(randomizers.poly.iter())
                .for_each(|(x, randomizer)| x.randomize_with(randomizer));
            Self::aggregate_randomized_signatures(&sig_randomized)
        }
    }
}

impl AddAssign for HotsSig {
    // Coefficient wise additions without mod reduction.
    fn add_assign(&mut self, other: Self) {
        // should not aggregate non-randomized signatures
        assert!(self.is_randomized);
        assert!(other.is_randomized);

        self.sigma
            .iter_mut()
            .zip(other.sigma.iter())
            .for_each(|(x, y)| *x = *x + *y)
    }
}
