mod ots_keys;
mod ots_param;
mod ots_sigs;

use crate::randomizer::Randomizers;
use crate::HOTSNTTPoly;
use crate::HOTSPoly;
use crate::HVCPoly;
use crate::Polynomial;
use crate::GAMMA;
use crate::PHI;
use crate::PHI_ALPHA_H;

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use sha2::Digest;

pub use ots_keys::HotsPK;
pub use ots_keys::RandomizedHOTSPK;
pub use ots_sigs::HotsSig;

pub struct HOTS;

// HOTS public parameters
#[derive(Debug, Clone, Copy)]
pub struct HotsParam {
    pub(crate) a: [HOTSNTTPoly; GAMMA],
}

// HOTS secret key
#[derive(Debug, Clone, Copy)]
pub struct HotsSK {
    pub(crate) s0: [HOTSNTTPoly; GAMMA],
    pub(crate) s1: [HOTSNTTPoly; GAMMA],
}

pub trait HomomorphicOneTimeSignature {
    type Param;
    type PK;
    type SK;
    type Signature;

    fn setup<R: Rng>(rng: &mut R) -> Self::Param;

    fn derive_sk(seed: &[u8; 32], counter: usize) -> Self::SK;

    fn key_gen(seed: &[u8; 32], counter: usize, pp: &Self::Param) -> (Self::PK, Self::SK);

    fn sign(sk: &Self::SK, message: &[u8]) -> Self::Signature;

    fn verify(pk: &Self::PK, message: &[u8], sig: &Self::Signature, pp: &Self::Param) -> bool;

    fn aggregate(sigs: &[Self::Signature], roots: &[HVCPoly]) -> Self::Signature;

    fn batch_verify(
        pks: &[Self::PK],
        message: &[u8],
        sig: &Self::Signature,
        roots: &[HVCPoly],
        pp: &Self::Param,
    ) -> bool;
}

impl HomomorphicOneTimeSignature for HOTS {
    type Param = HotsParam;
    type PK = HotsPK;
    type SK = HotsSK;
    type Signature = HotsSig;

    fn setup<R: Rng>(rng: &mut R) -> Self::Param {
        log::info!("HOTS key setup");
        let mut a = [HOTSNTTPoly::default(); GAMMA];
        a.iter_mut()
            .for_each(|x| *x = HOTSNTTPoly::from(&HOTSPoly::rand_poly(rng)));

        Self::Param { a }
    }

    fn derive_sk(seed: &[u8; 32], counter: usize) -> Self::SK {
        log::info!("HOTS derive sk");
        // initialize the rng with seed and counter
        let seed = [seed.as_ref(), counter.to_be_bytes().as_ref()].concat();
        let mut hasher = sha2::Sha256::new();
        hasher.update(seed);
        let seed = hasher.finalize().into();
        let mut rng = ChaCha20Rng::from_seed(seed);

        // sample the secret key
        let mut s0 = [HOTSNTTPoly::default(); GAMMA];
        let mut s1 = [HOTSNTTPoly::default(); GAMMA];

        s0.iter_mut()
            .for_each(|x| *x = HOTSNTTPoly::from(&HOTSPoly::rand_mod_p(&mut rng, PHI)));
        s1.iter_mut()
            .for_each(|x| *x = HOTSNTTPoly::from(&HOTSPoly::rand_mod_p(&mut rng, PHI_ALPHA_H)));
        Self::SK { s0, s1 }
    }

    fn key_gen(seed: &[u8; 32], counter: usize, pp: &Self::Param) -> (Self::PK, Self::SK) {
        log::info!("HOTS key generation");

        let sk = Self::derive_sk(seed, counter);

        // build the pk
        let mut pk = HotsPK::default();

        pp.a.iter()
            .zip(sk.s0.iter().zip(sk.s1.iter()))
            .for_each(|(&a, (&s0, &s1))| {
                pk.v0 += (&(a * s0)).into();
                pk.v1 += (&(a * s1)).into();
            });

        (pk, sk)
    }

    fn sign(sk: &Self::SK, message: &[u8]) -> Self::Signature {
        log::info!("HOTS signing");

        let mut sigma = [HOTSPoly::default(); GAMMA];
        let hm: HOTSNTTPoly = (&HOTSPoly::from_hash_message(message)).into();
        for (s, (&s0, &s1)) in sigma.iter_mut().zip(sk.s0.iter().zip(sk.s1.iter())) {
            *s = (&(s0 * hm + s1)).into();
        }
        HotsSig {
            sigma,
            is_randomized: false,
        }
    }

    fn verify(pk: &Self::PK, message: &[u8], sig: &Self::Signature, pp: &Self::Param) -> bool {
        log::info!("HOTS verification");

        //todo: check norm of signature
        let hm: HOTSNTTPoly = (&HOTSPoly::from_hash_message(message)).into();
        let mut left = HOTSNTTPoly::default();
        for (&a, s) in pp.a.iter().zip(sig.sigma.iter()) {
            left += a * HOTSNTTPoly::from(s)
        }
        let right = hm * HOTSNTTPoly::from(&pk.v0) + HOTSNTTPoly::from(&pk.v1);
        println!("left  {:?}\nright {:?}\n", left, right);
        // todo avoid the inverse NTT
        HOTSPoly::from(&left) == HOTSPoly::from(&right)
    }

    fn aggregate(sigs: &[Self::Signature], roots: &[HVCPoly]) -> Self::Signature {
        log::info!("HOTS aggregation");
        // check that length are correct
        assert_eq!(sigs.len(), roots.len());

        let randomizers = Randomizers::from_pks(roots);
        log::trace!("randomizers {:?}", randomizers);
        Self::Signature::aggregate_with_randomizers(sigs, &randomizers)
    }

    fn batch_verify(
        pks: &[Self::PK],
        message: &[u8],
        sig: &Self::Signature,
        roots: &[HVCPoly],
        pp: &Self::Param,
    ) -> bool {
        log::info!("HOTS batch verification");
        // check that length are correct
        assert_eq!(pks.len(), roots.len());

        let agg_pk = Self::PK::aggregate(pks, roots);
        Self::verify(&agg_pk, message, sig, pp)
    }
}

pub(crate) fn batch_verify_with_aggregated_pk(
    agg_pk: &RandomizedHOTSPK,
    message: &[u8],
    agg_sig: &HotsSig,
    pp: &HotsParam,
) -> bool {
    println!("agg hots pk in batch verify, before: {:?}", agg_pk);
    let agg_pk: HotsPK = agg_pk.into();
    println!("agg hots pk in batch verify, after: {:?}", agg_pk);
    HOTS::verify(&agg_pk, message, agg_sig, pp)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::poly::Polynomial;
    use rand::RngCore;

    #[test]
    fn test_hots() {
        let message = "this is the message to sign";
        let mut seed = [0u8; 32];
        let mut rng = ChaCha20Rng::from_seed(seed);
        let pp = HOTS::setup(&mut rng);

        for _ in 0..10 {
            rng.fill_bytes(&mut seed);
            let mut pks = Vec::new();
            let mut sigs = Vec::new();
            let mut roots = Vec::new();

            for counter in 0..100 {
                let (pk, sk) = HOTS::key_gen(&seed, counter, &pp);
                let sig = HOTS::sign(&sk, message.as_ref());
                assert!(HOTS::verify(&pk, message.as_ref(), &sig, &pp));
                pks.push(pk);
                sigs.push(sig);
                roots.push(HVCPoly::rand_poly(&mut rng));
            }
            let agg_sig = HOTS::aggregate(&sigs, &roots);
            assert!(HOTS::batch_verify(
                &pks,
                message.as_ref(),
                &agg_sig,
                &roots,
                &pp
            ));
            let pk_randomized: Vec<RandomizedHOTSPK> = pks.iter().map(|x| x.into()).collect();
            let agg_pk_randomized = RandomizedHOTSPK::aggregate(&pk_randomized, &roots);
            batch_verify_with_aggregated_pk(&agg_pk_randomized, message.as_ref(), &agg_sig, &pp);
        }
    }
}
