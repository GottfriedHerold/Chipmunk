use ark_std::{end_timer, start_timer};
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::poly::Polynomial;
use crate::{
    poly::HVCPoly, randomize_path::RandomizedPath, randomizer::Randomizers, HVCHash, HEIGHT,
};
use core::fmt;
use std::fmt::Display;

#[derive(Clone, Debug, PartialEq)]
pub struct Path {
    pub(crate) nodes: Vec<(HVCPoly, HVCPoly)>, // left and right nodes
    pub(crate) index: usize,
}

impl Default for Path {
    fn default() -> Self {
        Self {
            nodes: [(HVCPoly::default(), HVCPoly::default()); HEIGHT - 1].to_vec(),
            index: 0,
        }
    }
}

impl Display for Path {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let position_list = self.position_list();
        for ((i, (left, right)), is_right_node) in self.nodes.iter().enumerate().zip(position_list)
        {
            writeln!(
                f,
                "{}-th node: \nleft    {}\nright   {}\nis right {}",
                i, left, right, is_right_node
            )?;
        }

        Ok(())
    }
}

impl Path {
    /// The position of on_path node in `leaf_and_sibling_hash` and `non_leaf_and_sibling_hash_path`.
    /// `position[i]` is 0 (false) iff `i`th on-path node from top to bottom is on the left.
    ///
    /// This function simply converts `self.leaf_index` to boolean array in big endian form.
    fn position_list(&'_ self) -> impl '_ + Iterator<Item = bool> {
        (0..self.nodes.len() + 1)
            .map(move |i| ((self.index >> i) & 1) != 0)
            .rev()
    }

    /// verifies the path against a root
    pub fn verify(&self, root: &HVCPoly, hasher: &HVCHash) -> bool {
        let timer = start_timer!(|| "hvc path verify");
        // check that the first two elements hashes to root
        if hasher.decom_then_hash(&self.nodes[0].0, &self.nodes[0].1) != *root {
            log::error!("hasher does not matcht the root");
            return false;
        }
        let position_list: Vec<_> = self.position_list().collect();

        let checks = self
            .nodes
            .clone()
            .into_par_iter()
            .enumerate()
            .skip(1)
            .map(|(i, (left, right))| {
                if position_list[i] {
                    hasher.decom_then_hash(&left, &right) == self.nodes[i - 1].1
                } else {
                    hasher.decom_then_hash(&left, &right) == self.nodes[i - 1].0
                }
            })
            .collect::<Vec<_>>()
            .iter()
            .fold(true, |acc, mk| acc && *mk);

        end_timer!(timer);
        checks
    }

    pub(crate) fn aggregate_with_randomizers(
        paths: &[Self],
        randomizers: &Randomizers,
    ) -> RandomizedPath {
        log::info!("hvc aggregation with randomizers");
        // check that length are correct
        let len = paths.len();
        assert_eq!(len, randomizers.poly.len());
        // check that we aggregate for a same index
        for e in paths.iter().skip(1) {
            assert_eq!(e.index, paths[0].index)
        }

        let randomized_paths: Vec<RandomizedPath> = paths.iter().map(|x| x.into()).collect();
        RandomizedPath::aggregate_with_randomizers(&randomized_paths, randomizers)
    }

    /// Aggregate a set of paths
    pub fn aggregation(paths: &[Self], roots: &[HVCPoly]) -> RandomizedPath {
        log::info!("hvc aggregation");
        // get and apply the randomizers
        let randomizers = Randomizers::from_pks(roots);
        Self::aggregate_with_randomizers(paths, &randomizers)
    }

    pub fn random_for_testing<R: rand::Rng>(rng: &mut R, hasher: &HVCHash) -> (Self, HVCPoly) {
        let mut nodes = vec![(HVCPoly::rand_poly(rng), HVCPoly::rand_poly(rng))];

        for i in 1..HEIGHT - 1 {
            let left = hasher.decom_then_hash(&nodes[i - 1].0, &nodes[i - 1].1);
            nodes.push((left, HVCPoly::rand_poly(rng)))
        }
        let root = hasher.decom_then_hash(&nodes[HEIGHT - 2].0, &nodes[HEIGHT - 2].1);
        nodes.reverse();
        let mut path = Self::default();
        path.nodes.clone_from_slice(&nodes);
        path.index = 0;

        (path, root)
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_path() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let hasher = HVCHash::init(&mut rng);
        for _ in 0..100 {
            let (path, root) = Path::random_for_testing(&mut rng, &hasher);
            assert!(path.verify(&root, &hasher))
        }
    }

    #[test]
    fn test_path_conversion() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let hasher = HVCHash::init(&mut rng);
        for _ in 0..100 {
            let (path, _root) = Path::random_for_testing(&mut rng, &hasher);
            let randomize_path: RandomizedPath = (&path).into();
            let path_rec = (&randomize_path).into();
            assert_eq!(path, path_rec)
        }
    }
}
