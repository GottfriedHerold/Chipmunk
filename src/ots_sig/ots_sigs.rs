use ark_std::{end_timer, start_timer};
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::{
    io::{Read, Write},
    ops::AddAssign,
};

use crate::{poly::HOTSPoly, randomizer::Randomizers, TerPolyCoeffEncoding, GAMMA};

// HOTS signature
#[derive(Debug, Clone, Copy, PartialEq, Default)]
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
            let timer = start_timer!(|| format!("HOTS aggregate {} signatures", sigs.len()));

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
            let res = Self::aggregate_randomized_signatures(&sig_randomized);
            end_timer!(timer);
            res
        }
        #[cfg(not(feature = "parallel"))]
        {
            let timer = start_timer!(|| format!("HOTS aggregate {} signatures", sigs.len()));

            let mut sig_randomized = sigs.to_vec();
            sig_randomized
                .iter_mut()
                .zip(randomizers.poly.iter())
                .for_each(|(x, randomizer)| x.randomize_with(randomizer));
            let res = Self::aggregate_randomized_signatures(&sig_randomized);
            end_timer!(timer);
            res
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

// (De)Serializations
impl HotsSig {
    pub(crate) fn serialize<W: Write>(&self, mut writer: W, is_aggregated: bool) {
        let timer = start_timer!(|| "HotsSig serialization");
        assert_eq!(is_aggregated, self.is_randomized);

        self.sigma.iter().for_each(|poly| {
            poly.coeffs
                .iter()
                .for_each(|coeff| writer.write_all(coeff.to_le_bytes().as_ref()).unwrap())
        });

        writer
            .write_all((self.is_randomized as u8).to_le_bytes().as_ref())
            .unwrap();

        end_timer!(timer);
    }

    pub(crate) fn deserialize<R: Read>(mut reader: R) -> Self {
        let timer = start_timer!(|| "HotsSig deserialization");
        let mut res = Self::default();
        let mut buf4 = [0u8; 4];
        let mut buf = [0u8];

        res.sigma.iter_mut().for_each(|poly| {
            poly.coeffs.iter_mut().for_each(|coeff| {
                reader.read_exact(&mut buf4).unwrap();
                *coeff = i32::from_le_bytes(buf4);
            })
        });

        reader.read_exact(&mut buf).unwrap();
        assert!(buf[0] == 0 || buf[0] == 1);
        res.is_randomized = buf[0] != 0;

        end_timer!(timer);
        res
    }
}
