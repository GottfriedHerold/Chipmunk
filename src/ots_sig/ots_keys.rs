use std::{
    io::{Read, Write},
    ops::AddAssign,
};

use crate::{
    poly::{HOTSPoly, TerPolyCoeffEncoding},
    randomizer::Randomizers,
    HOTSHash, HVCPoly, HOTS_WIDTH,
};
use ark_std::{end_timer, start_timer};
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefIterator};

// HOTS public key
#[derive(Debug, Default, Clone, Copy)]
pub struct HotsPK {
    pub(crate) v0: HOTSPoly,
    pub(crate) v1: HOTSPoly,
}

impl HotsPK {
    pub(crate) fn digest(&self, hasher: &HOTSHash) -> HVCPoly {
        hasher.hash_separate_inputs(&self.v0.decompose_r(), &self.v1.decompose_r())
    }

    /// Aggregate multiple PKs into a single PK
    pub(crate) fn aggregate(pks: &[Self], roots: &[HVCPoly]) -> Self {
        // get and apply the randomizers
        let randomizers = Randomizers::from_pks(roots);
        Self::aggregate_with_randomizers(pks, &randomizers)
    }

    /// Aggregate a set of pks with randomizes
    pub(crate) fn aggregate_with_randomizers(pks: &[Self], randomizers: &Randomizers) -> Self {
        let mut pk_and_randomizer: Vec<(Self, HOTSPoly)> = pks
            .iter()
            .zip(randomizers.poly.iter())
            .map(|(&pk, r)| (pk, HOTSPoly::from(r)))
            .collect();

        #[cfg(feature = "parallel")]
        pk_and_randomizer.par_iter_mut().for_each(|(pk, r)| {
            pk.v0 = pk.v0 * *r;
            pk.v1 = pk.v1 * *r;
        });

        #[cfg(not(feature = "parallel"))]
        pk_and_randomizer.iter_mut().for_each(|(pk, r)| {
            pk.v0 = pk.v0 * *r;
            pk.v1 = pk.v1 * *r;
        });

        let mut agg_pk = pk_and_randomizer[0].0;

        for (pk, _r) in pk_and_randomizer.iter().skip(1) {
            agg_pk.v0 += pk.v0;
            agg_pk.v1 += pk.v1;
        }
        agg_pk
    }
}

// HOTS public key
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct RandomizedHOTSPK {
    pub(crate) v0: [HVCPoly; HOTS_WIDTH],
    pub(crate) v1: [HVCPoly; HOTS_WIDTH],
    pub(crate) is_randomized: bool,
}

impl From<&HotsPK> for RandomizedHOTSPK {
    fn from(pk: &HotsPK) -> Self {
        log::info!("randomize hots pk");
        RandomizedHOTSPK {
            v0: pk.v0.decompose_r(),
            v1: pk.v1.decompose_r(),
            is_randomized: false,
        }
    }
}

impl From<&RandomizedHOTSPK> for HotsPK {
    fn from(pk: &RandomizedHOTSPK) -> Self {
        HotsPK {
            v0: HOTSPoly::projection_r(&pk.v0),
            v1: HOTSPoly::projection_r(&pk.v1),
        }
    }
}

impl RandomizedHOTSPK {
    pub(crate) fn randomize_with(&mut self, ternary_coeffs: &TerPolyCoeffEncoding) {
        if self.is_randomized {
            panic!("already randomized")
        }

        #[cfg(not(feature = "parallel"))]
        self.v0.iter_mut().for_each(|x| {
            *x = SignedPoly::ter_mul_bin(&ternary_coeffs, x);
        });
        #[cfg(not(feature = "parallel"))]
        self.v1.iter_mut().for_each(|x| {
            *x = SignedPoly::ter_mul_bin(&ternary_coeffs, x);
        });

        #[cfg(feature = "parallel")]
        self.v0.par_iter_mut().for_each(|x| {
            *x = HVCPoly::ter_mul_bin(&ternary_coeffs, x);
        });
        #[cfg(feature = "parallel")]
        self.v1.par_iter_mut().for_each(|x| {
            *x = HVCPoly::ter_mul_bin(&ternary_coeffs, x);
        });

        self.is_randomized = true;
    }

    pub(crate) fn digest(&self, hasher: &HOTSHash) -> HVCPoly {
        hasher.hash_separate_inputs(&self.v0, &self.v1)
    }

    /// Aggregate multiple PKs into a single PK
    pub(crate) fn aggregate(pks: &[Self], roots: &[HVCPoly]) -> Self {
        // get and apply the randomizers
        let randomizers = Randomizers::from_pks(roots);
        Self::aggregate_with_randomizers(pks, &randomizers)
    }

    /// Aggregate a set of pks with randomizes
    pub(crate) fn aggregate_with_randomizers(pks: &[Self], randomizers: &Randomizers) -> Self {
        let timer = start_timer!(|| format!("aggregating {} HOTS pks", pks.len()));
        let mut randomized_pks: Vec<Self> = pks.to_vec();
        randomized_pks
            .par_iter_mut()
            .zip(randomizers.poly.par_iter())
            .for_each(|(x, r)| x.randomize_with(r));

        let mut res = randomized_pks[0];
        randomized_pks.iter().skip(1).for_each(|x| res += *x);
        end_timer!(timer);
        res
    }
}

impl AddAssign for RandomizedHOTSPK {
    // Coefficient wise additions without mod reduction.
    fn add_assign(&mut self, other: Self) {
        self.v0
            .iter_mut()
            .zip(other.v0.iter())
            .for_each(|(x, y)| *x += *y);
        self.v1
            .iter_mut()
            .zip(other.v1.iter())
            .for_each(|(x, y)| *x += *y);
    }
}

// (De)Serializations
impl RandomizedHOTSPK {
    pub(crate) fn serialize<W: Write>(&self, mut writer: W, is_aggregated: bool) {
        let timer = start_timer!(|| "RandomizedHOTSPK serialization");
        assert_eq!(is_aggregated, self.is_randomized);

        self.v0.iter().for_each(|poly| {
            poly.coeffs
                .iter()
                .for_each(|coeff| writer.write_all(coeff.to_le_bytes().as_ref()).unwrap())
        });
        self.v1.iter().for_each(|poly| {
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
        let timer = start_timer!(|| "RandomizedHOTSPK deserialization");
        let mut res = Self::default();
        let mut buf4 = [0u8; 4];
        let mut buf = [0u8];

        res.v0.iter_mut().for_each(|poly| {
            poly.coeffs.iter_mut().for_each(|coeff| {
                reader.read_exact(&mut buf4).unwrap();
                *coeff = i32::from_le_bytes(buf4);
            })
        });
        res.v1.iter_mut().for_each(|poly| {
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

#[cfg(test)]
mod tests {

    use std::io::Cursor;

    use super::*;
    use crate::poly::Polynomial;
    use crate::{HomomorphicOneTimeSignature, HOTS};
    use rand::{RngCore, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_homomorphic_hash() {
        let mut seed = [0u8; 32];
        let mut rng = ChaCha20Rng::from_seed(seed);
        let pp = HOTS::setup(&mut rng);
        let hasher = HOTSHash::init(&mut rng);

        for _ in 0..10 {
            rng.fill_bytes(&mut seed);
            let mut pks = Vec::new();
            let mut pks_randomized = Vec::new();
            let mut roots = Vec::new();
            let mut digests = Vec::new();

            for counter in 0..1 {
                let (pk, _sk) = HOTS::key_gen(&seed, counter, &pp);
                let rand_pk = RandomizedHOTSPK::from(&pk);
                let digest = rand_pk.digest(&hasher);

                pks_randomized.push(rand_pk);
                pks.push(pk);
                digests.push(digest);
                roots.push(HVCPoly::rand_poly(&mut rng));
            }
            let randomizers = Randomizers::from_pks(&roots);
            let agg_pk_randomized = RandomizedHOTSPK::aggregate(&pks_randomized, &roots);
            {
                let mut bytes = vec![];
                agg_pk_randomized.serialize(&mut bytes, true);
                let buf = Cursor::new(bytes);
                let agg_pk_reconstructed = RandomizedHOTSPK::deserialize(buf);
                assert_eq!(agg_pk_randomized, agg_pk_reconstructed);
            }

            let agg_digest = agg_pk_randomized.digest(&hasher);
            let mut agg_digest_rec = HVCPoly::default();

            for (&d, r) in digests.iter().zip(randomizers.poly.iter()) {
                agg_digest_rec += d * HVCPoly::from(r);
            }
            assert_eq!(agg_digest, agg_digest_rec);
        }
    }
}
