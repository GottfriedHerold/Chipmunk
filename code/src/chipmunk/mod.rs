use crate::path::Path;
use crate::randomize_path::RandomizedPath;
use crate::randomizer::Randomizers;
use crate::Tree;
use crate::{
    batch_verify_with_aggregated_pk, HOTSHash, HVCHash, HVCPoly, HomomorphicOneTimeSignature,
    HotsParam, HotsSig, MultiSig, RandomizedHOTSPK, HEIGHT, HOTS,
};
use ark_std::{end_timer, start_timer};
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

#[cfg(test)]
mod tests;

pub struct Chipmunk;

#[derive(Debug, Clone)]
pub struct ChipmunkParam {
    hvc_hasher: HVCHash,
    hots_hasher: HOTSHash,
    hots_param: HotsParam,
}
#[derive(Debug, Clone)]
pub struct ChipmunkSK {
    sk_seed: [u8; 32],
    tree: Tree,
}

pub type ChipmunkPK = HVCPoly;
#[derive(Debug, Clone)]
pub struct ChipmunkSignature {
    path: RandomizedPath,
    hots_pk: RandomizedHOTSPK,
    hots_sig: HotsSig,
}

impl MultiSig for Chipmunk {
    type Param = ChipmunkParam;
    type PK = ChipmunkPK;
    type SK = ChipmunkSK;
    type Signature = ChipmunkSignature;

    fn setup<R: Rng>(rng: &mut R) -> Self::Param {
        Self::Param {
            hvc_hasher: HVCHash::init(rng),
            hots_hasher: HOTSHash::init(rng),
            hots_param: HOTS::setup(rng),
        }
    }

    fn key_gen(seed: &[u8; 32], pp: &Self::Param) -> (Self::PK, Self::SK) {
        let timer = start_timer!(|| "Chipmunk key gen");

        let leaf_timer = start_timer!(|| "leaf generation");
        let mut pk_digests = vec![HVCPoly::default(); 1 << (HEIGHT - 1)];

        #[cfg(not(feature = "parallel"))]
        pk_digests.iter_mut().enumerate().for_each(|(index, pkd)| {
            let (pk, _sk) = HOTS::key_gen(seed, index, &pp.hots_param);
            *pkd = pk.digest(&pp.hots_hasher)
        });

        #[cfg(feature = "parallel")]
        pk_digests
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, pkd)| {
                let (pk, _sk) = HOTS::key_gen(seed, index, &pp.hots_param);
                *pkd = pk.digest(&pp.hots_hasher)
            });
        end_timer!(leaf_timer);

        let tree = Tree::new_with_leaf_nodes(&pk_digests, &pp.hvc_hasher);
        end_timer!(timer);
        (
            tree.root(),
            ChipmunkSK {
                sk_seed: *seed,
                tree,
            },
        )
    }

    fn sign(sk: &Self::SK, index: usize, message: &[u8], pp: &Self::Param) -> Self::Signature {
        let timer = start_timer!(|| "Chipmunk Signing");
        let path = sk.tree.gen_proof(index);
        let (hots_pk, hots_sk) = HOTS::key_gen(&sk.sk_seed, index, &pp.hots_param);

        let hots_sig = HOTS::sign(&hots_sk, message);
        let res = ChipmunkSignature {
            path: (&path).into(),
            hots_pk: (&hots_pk).into(),
            hots_sig,
        };
        end_timer!(timer);
        res
    }

    fn verify(pk: &Self::PK, message: &[u8], sig: &Self::Signature, pp: &Self::Param) -> bool {
        let timer = start_timer!(|| "Chipmunk verify");
        // check signature against hots pk
        let hots_pk = (&sig.hots_pk).into();

        if !HOTS::verify(&hots_pk, message, &sig.hots_sig, &pp.hots_param) {
            log::error!("hots verification failed");
            return false;
        }

        // check hots public key membership
        let path = Path::from(&sig.path);
        if !path.verify(pk, &pp.hvc_hasher) {
            log::error!("hvc verification failed");
            return false;
        }

        let pk_digest = hots_pk.digest(&pp.hots_hasher);

        let res = if sig.path.index & 1 == 0 {
            pk_digest == path.nodes[HEIGHT - 2].0
        } else {
            pk_digest == path.nodes[HEIGHT - 2].1
        };
        if !res {
            log::error!("pk does not match the node");
        }
        end_timer!(timer);
        res
    }

    fn aggregate(sigs: &[Self::Signature], roots: &[HVCPoly]) -> Self::Signature {
        let timer = start_timer!(|| format!("aggregating {} signatures", sigs.len()));
        let randomizers = Randomizers::from_pks(roots);

        // aggregate HOTS pk
        let pks: Vec<RandomizedHOTSPK> = sigs.iter().map(|x| x.hots_pk).collect();
        let agg_pk = RandomizedHOTSPK::aggregate_with_randomizers(&pks, &randomizers);

        // aggregate HOTS sig
        let hots_sigs: Vec<HotsSig> = sigs.iter().map(|x| x.hots_sig).collect();
        let agg_sig = HotsSig::aggregate_with_randomizers(&hots_sigs, &randomizers);

        // aggregate the membership proof
        let membership_proofs: Vec<RandomizedPath> = sigs.iter().map(|x| x.path.clone()).collect();
        let agg_proof =
            RandomizedPath::aggregate_with_randomizers(&membership_proofs, &randomizers);

        end_timer!(timer);
        Self::Signature {
            path: agg_proof,
            hots_pk: agg_pk,
            hots_sig: agg_sig,
        }
    }

    fn batch_verify(
        pks: &[Self::PK],
        message: &[u8],
        sig: &Self::Signature,
        pp: &Self::Param,
    ) -> bool {
        let timer = start_timer!(|| format!("Chipmunk batch verify {} signatures", pks.len()));
        if !batch_verify_with_aggregated_pk(&sig.hots_pk, message, &sig.hots_sig, &pp.hots_param) {
            log::error!("HOTS batch verification failed");
            return false;
        }
        if !sig.path.verify(pks, &pp.hvc_hasher) {
            log::error!("Path batch verification failed");
            return false;
        }
        let res = if sig.path.index & 1 == 0 {
            sig.hots_pk.digest(&pp.hots_hasher)
                == HVCPoly::projection(&sig.path.nodes[HEIGHT - 2].0)
        } else {
            sig.hots_pk.digest(&pp.hots_hasher)
                == HVCPoly::projection(&sig.path.nodes[HEIGHT - 2].1)
        };
        end_timer!(timer);
        res
    }
}
