// =================SHARED PARAM===============
/// degree of polynomial
pub const N: usize = 512;
/// Security parameter
pub const SEC_PARAM: usize = 112;
/// non-zero entries in randomizer
pub const ALPHA: usize = 16;
/// Height of the tree
pub const HEIGHT: usize = 5;
/// base of the decomposition: the decomposed poly have coefficients between [-zeta, zeta]
pub const ZETA: u32 = 29;
/// the arity: 2 * zeta + 1
pub const TWO_ZETA_PLUS_ONE: u32 = 59;

// =================HOTS PARAM=================
/// q for large ring, HOTS modulus 3168257
pub const HOTS_MODULUS: i32 = 3168257;
/// 1/N mod q
pub const HOTS_ONE_OVER_N: i32 = 3162069;
/// (q-1)/2 = 1584128
pub const HOTS_MODULUS_OVER_TWO: i32 = 1584128;
/// number of ring elements during decomposition
/// require (2 * zeta + 1)^HOTS_WIDTH > HOTS_MODULUS
pub const HOTS_WIDTH: usize = 4;
/// the largest multiple of q that is smaller than 2^32
pub const HOTS_SAMPLE_THRESHOLD: u32 = 4292988235;
/// the number of polynomials in a decomposed poly
pub const GAMMA: usize = 6;
/// hamming weight of the hash of the message
pub const ALPHA_H: usize = 37;
/// the infinity norm bound for s0
/// i.e., a fresh s_0 has coefficients in [-PHI, PHI]
pub const PHI: u32 = 13;
///the largest multiple of (2 * phi + 1) that is smaller than 2^32
pub const PHI_SAMPLE_THRESHOLD: u32 = 4294967274;
/// norm bound of s_1 = phi * alpha_H
/// i.e., a fresh s_1 has coefficients in [-phi * alpha_H, phi * alpha_H]
pub const PHI_ALPHA_H: u32 = 481;
///the largest multiple of (2 * PHI_ALPHA_H + 1) that is smaller than 2^32
pub const PHI_ALPHA_H_SAMPLE_THRESHOLD: u32 = 4294966518;

// =================HVC PARAM==================
/// q for small ring, HVC modulus 202753
pub const HVC_MODULUS: i32 = 202753;
/// 1/N mod q
pub const HVC_ONE_OVER_N: i32 = 202357;
/// (q-1)/2
pub const HVC_MODULUS_OVER_TWO: i32 = 101376;
/// the largest multiple of q that is smaller than 2^32
pub const HVC_SAMPLE_THRESHOLD: u32 = 4294916799;
/// number of ring elements during decomposition
/// require (2 * zeta + 1)^HVC_WIDTH > HVC_MODULUS
pub const HVC_WIDTH: usize = 3;
