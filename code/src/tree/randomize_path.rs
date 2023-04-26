use ark_std::{end_timer, start_timer};
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::path::Path;
use crate::{
    poly::{HVCPoly, TerPolyCoeffEncoding},
    randomizer::Randomizers,
    HVCHash, HEIGHT, HVC_WIDTH,
};
use core::fmt;
use std::fmt::Display;
use std::ops::{Add, AddAssign};

#[derive(Clone, Debug)]
pub struct RandomizedPath {
    pub(crate) nodes: Vec<([HVCPoly; HVC_WIDTH], [HVCPoly; HVC_WIDTH])>,
    pub(crate) index: usize,
    pub(crate) is_randomized: bool,
}

impl Default for RandomizedPath {
    fn default() -> Self {
        Self {
            nodes: [(
                [HVCPoly::default(); HVC_WIDTH],
                [HVCPoly::default(); HVC_WIDTH],
            ); HEIGHT - 1]
                .to_vec(),
            index: 0,
            is_randomized: false,
        }
    }
}

impl Display for RandomizedPath {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let position_list = self.position_list();
        for ((i, (left, right)), is_right_node) in self.nodes.iter().enumerate().zip(position_list)
        {
            writeln!(f, "{}-th node: is right {}", i, is_right_node)?;
            for (l, r) in left.iter().zip(right.iter()) {
                writeln!(f, "{}", l)?;
                writeln!(f, "{}", r)?;
            }
        }

        Ok(())
    }
}

// todo: improve the code to avoid `clone_from_slice`
impl From<&Path> for RandomizedPath {
    fn from(p: &Path) -> Self {
        // seems that the overhead of parallelizing the conversion is enormous
        // and we are better off without parallelization.
        //
        // #[cfg(feature = "parallel")]
        // let nodes: Vec<(Vec<HVCPoly>, Vec<HVCPoly>)> = p
        //     .nodes
        //     .clone()
        //     .into_par_iter()
        //     .map(|(left, right)| (left.decompose(), right.decompose()))
        //     .collect();
        //
        // #[cfg(not(feature = "parallel"))]
        let nodes: Vec<_> = p
            .nodes
            .iter()
            .map(|(left, right)| (left.decompose(), right.decompose()))
            .collect();
        let mut res = Self::default();
        res.nodes.clone_from_slice(&nodes);
        res.index = p.index;
        res.is_randomized = false;

        res
    }
}

impl From<&RandomizedPath> for Path {
    fn from(r: &RandomizedPath) -> Self {
        // #[cfg(feature = "parallel")]
        // let nodes = r
        //     .nodes
        //     .clone()
        //     .into_par_iter()
        //     .map(|(left, right)| (HVCPoly::projection(&left), HVCPoly::projection(&right)))
        //     .collect();
        //
        // #[cfg(not(feature = "parallel"))]
        let nodes: Vec<_> = r
            .nodes
            .iter()
            .map(|(left, right)| (HVCPoly::projection(left), HVCPoly::projection(right)))
            .collect();
        let mut res = Self::default();
        res.nodes.clone_from_slice(&nodes);
        res.index = r.index;
        res
    }
}

impl RandomizedPath {
    /// The position of on_path node in `leaf_and_sibling_hash` and `non_leaf_and_sibling_hash_path`.
    /// `position[i]` is 0 (false) iff `i`th on-path node from top to bottom is on the left.
    ///
    /// This function simply converts `self.leaf_index` to boolean array in big endian form.
    fn position_list(&'_ self) -> impl '_ + Iterator<Item = bool> {
        (0..self.nodes.len() + 1)
            .map(move |i| ((self.index >> i) & 1) != 0)
            .rev()
    }

    pub fn randomize_with(&mut self, ternary_coeffs: &TerPolyCoeffEncoding) {
        if self.is_randomized {
            panic!("already randomized")
        }

        #[cfg(not(feature = "parallel"))]
        self.nodes.iter_mut().for_each(|(left, right)| {
            left.iter_mut()
                .for_each(|x| *x = HVCPoly::ter_mul_bin(&ternary_coeffs, x));
            right
                .iter_mut()
                .for_each(|x| *x = HVCPoly::ter_mul_bin(&ternary_coeffs, x));
        });

        #[cfg(feature = "parallel")]
        self.nodes.par_iter_mut().for_each(|(left, right)| {
            left.par_iter_mut()
                .for_each(|x| *x = HVCPoly::ter_mul_bin(&ternary_coeffs, x));
            right
                .par_iter_mut()
                .for_each(|x| *x = HVCPoly::ter_mul_bin(&ternary_coeffs, x));
        });

        self.is_randomized = true;
    }

    pub(crate) fn aggregate_with_randomizers(paths: &[Self], randomizers: &Randomizers) -> Self {
        let timer = start_timer!(|| format!("aggregate {} paths", paths.len()));
        let mut randomized_paths: Vec<RandomizedPath> = paths.to_vec();
        randomized_paths
            .iter_mut()
            .zip(randomizers.poly.iter())
            .for_each(|(path, randomizer)| path.randomize_with(randomizer));

        // aggregate the result
        let mut res = randomized_paths[0].clone();

        randomized_paths
            .iter()
            .skip(1)
            .for_each(|target| res = &res + target);
        end_timer!(timer);
        res
    }

    /// verifies the path against a list of root
    pub fn verify(&self, roots: &[HVCPoly], hasher: &HVCHash) -> bool {
        let timer = start_timer!(|| format!("batch verify {} paths", roots.len()));
        // recompute the root
        let randomziers = Randomizers::from_pks(roots);
        let mut root = HVCPoly::default();
        roots
            .iter()
            .zip(randomziers.poly.iter())
            .for_each(|(&rt, rand)| root += HVCPoly::from(rand) * rt);

        // check that the first two elements hashes to root
        if hasher.hash_separate_inputs(self.nodes[0].0.as_ref(), self.nodes[0].1.as_ref()) != root {
            return false;
        }
        let position_list = self.position_list();

        for ((i, (left, right)), is_right_node) in
            self.nodes.iter().enumerate().zip(position_list).skip(1)
        {
            let digest = hasher.hash_separate_inputs(left, right);
            if is_right_node {
                if digest != HVCPoly::projection(&self.nodes[i - 1].1) {
                    return false;
                }
            } else if digest != HVCPoly::projection(&self.nodes[i - 1].0) {
                return false;
            }
        }
        end_timer!(timer);
        true
    }
}

impl Add for &RandomizedPath {
    type Output = RandomizedPath;

    // Coefficient wise additions without mod reduction.
    fn add(self, other: Self) -> RandomizedPath {
        assert_eq!(self.index, other.index);

        let mut res = self.clone();
        res.nodes
            .iter_mut()
            .zip(other.nodes.iter())
            .for_each(|(x, y)| {
                x.0.iter_mut()
                    .zip(y.0.iter())
                    .for_each(|(x0, y0)| *x0 += *y0);
                x.1.iter_mut()
                    .zip(y.1.iter())
                    .for_each(|(x1, y1)| *x1 += *y1);
            });

        res
    }
}

impl AddAssign for RandomizedPath {
    // Coefficient wise additions without mod reduction.
    fn add_assign(&mut self, other: Self) {
        assert_eq!(self.index, other.index);

        self.nodes
            .iter_mut()
            .zip(other.nodes.iter())
            .for_each(|(x, y)| {
                x.0.iter_mut()
                    .zip(y.0.iter())
                    .for_each(|(x0, y0)| *x0 += *y0);
                x.1.iter_mut()
                    .zip(y.1.iter())
                    .for_each(|(x1, y1)| *x1 += *y1);
            });
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_randomized_path() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        for _ in 0..10 {
            let hasher = HVCHash::init(&mut rng);
            let mut paths = vec![];
            let mut roots = vec![];
            for _ in 0..100 {
                let (path, root) = Path::random_for_testing(&mut rng, &hasher);
                assert!(path.verify(&root, &hasher));
                paths.push(path);
                roots.push(root);
            }

            let path = Path::aggregation(&paths, &roots);
            assert!(path.verify(&roots, &hasher))
        }
    }
}
