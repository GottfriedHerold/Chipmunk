mod chipmunk;
mod hash;
mod ots_sig;
mod param;
mod poly;
mod randomizer;
mod tree;

pub use chipmunk::*;
pub use hash::*;
pub use ots_sig::*;
pub use param::*;
pub use poly::*;
use rand::Rng;
pub use tree::*;

pub trait MultiSig {
    type Param;
    type PK;
    type SK;
    type FreshSignature;
    type AggregatedSignature;

    fn setup<R: Rng>(rng: &mut R) -> Self::Param;

    fn key_gen(seed: &[u8; 32], pp: &Self::Param) -> (Self::PK, Self::SK);

    /// Sign a message for the `index` time slot
    fn sign(sk: &Self::SK, index: usize, message: &[u8], pp: &Self::Param) -> Self::FreshSignature;

    fn verify(pk: &Self::PK, message: &[u8], sig: &Self::FreshSignature, pp: &Self::Param) -> bool;

    fn aggregate(sigs: &[Self::FreshSignature], roots: &[HVCPoly]) -> Self::AggregatedSignature;

    fn batch_verify(
        pks: &[Self::PK],
        message: &[u8],
        sig: &Self::AggregatedSignature,
        pp: &Self::Param,
    ) -> bool;
}
