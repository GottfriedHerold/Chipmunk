use ark_std::{end_timer, start_timer};
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator};

use crate::path::Path;
use crate::{
    poly::{HVCPoly, TerPolyCoeffEncoding},
    randomizer::Randomizers,
    HVCHash, HEIGHT, HVC_WIDTH,
};
use core::fmt;
use std::fmt::Display;
use std::io::{Read, Write};
use std::ops::{Add, AddAssign};

#[derive(Clone, Debug, PartialEq)]
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
            .map(|(left, right)| (left.decompose_r(), right.decompose_r()))
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
            .map(|(left, right)| (HVCPoly::projection_r(left), HVCPoly::projection_r(right)))
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
            .par_iter_mut()
            .zip(randomizers.poly.par_iter())
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
    pub fn complete(&mut self, roots: &[HVCPoly], hasher: &HVCHash) {
        let timer = start_timer!(|| format!("complete {} paths", roots.len()));

        // recompute the root
        let randomziers = Randomizers::from_pks(roots);

        let tmp: Vec<_> = roots
            .par_iter()
            .zip(randomziers.poly.par_iter())
            .map(|(&rt, rand)| HVCPoly::from(rand) * rt)
            .collect();
        let root = tmp.iter().fold(HVCPoly::default(), |acc, mk| acc + *mk);

        // check that the first two elements hashes to root
        self.nodes[0].1[HVC_WIDTH - 1] =
            hasher.derive_missing_input(self.nodes[0].0.as_ref(), self.nodes[0].1.as_ref(), &root);

        let position_list: Vec<_> = self.position_list().collect();

        // todo: parallelization?
        for i in 1..self.nodes.len() {
            if position_list[i] {
                self.nodes[i].1[HVC_WIDTH - 1] = hasher.derive_missing_input(
                    self.nodes[i].0.as_ref(),
                    self.nodes[i].1.as_ref(),
                    &HVCPoly::projection_r(&self.nodes[i - 1].1),
                );
            } else {
                self.nodes[i].1[HVC_WIDTH - 1] = hasher.derive_missing_input(
                    self.nodes[i].0.as_ref(),
                    self.nodes[i].1.as_ref(),
                    &HVCPoly::projection_r(&self.nodes[i - 1].0),
                );
            }
        }
        end_timer!(timer);
    }

    /// verifies the path against a list of root
    pub fn verify(&self, roots: &[HVCPoly], hasher: &HVCHash) -> bool {
        let timer = start_timer!(|| format!("batch verify {} paths", roots.len()));
        // recompute the root
        let randomziers = Randomizers::from_pks(roots);

        let tmp: Vec<_> = roots
            .par_iter()
            .zip(randomziers.poly.par_iter())
            .map(|(&rt, rand)| HVCPoly::from(rand) * rt)
            .collect();
        let root = tmp.iter().fold(HVCPoly::default(), |acc, mk| acc + *mk);

        // check that the first two elements hashes to root
        if hasher.hash_separate_inputs(self.nodes[0].0.as_ref(), self.nodes[0].1.as_ref()) != root {
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
                    hasher.hash_separate_inputs(&left, &right)
                        == HVCPoly::projection_r(&self.nodes[i - 1].1)
                } else {
                    hasher.hash_separate_inputs(&left, &right)
                        == HVCPoly::projection_r(&self.nodes[i - 1].0)
                }
            })
            .collect::<Vec<_>>()
            .iter()
            .fold(true, |acc, mk| acc && *mk);

        end_timer!(timer);
        checks
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

// (De)Serializations
impl RandomizedPath {
    pub(crate) fn serialize<W: Write>(
        &self,
        mut writer: W,
        is_aggregated: bool,
        skip_last_poly: bool,
    ) {
        let timer = start_timer!(|| "RandomizedPath serialization");
        assert_eq!(is_aggregated, self.is_randomized);

        self.nodes.iter().for_each(|(left, right)| {
            left.iter().for_each(|poly| {
                poly.coeffs
                    .iter()
                    .for_each(|coeff| writer.write_all(coeff.to_le_bytes().as_ref()).unwrap())
            });
            right
                .iter()
                .take(HVC_WIDTH - skip_last_poly as usize)
                .for_each(|poly| {
                    poly.coeffs
                        .iter()
                        .for_each(|coeff| writer.write_all(coeff.to_le_bytes().as_ref()).unwrap())
                });
        });
        writer.write_all(self.index.to_le_bytes().as_ref()).unwrap();
        writer
            .write_all((self.is_randomized as u8).to_le_bytes().as_ref())
            .unwrap();

        end_timer!(timer);
    }

    pub(crate) fn deserialize<R: Read>(mut reader: R, skip_last_poly: bool) -> Self {
        let timer = start_timer!(|| "RandomizedPath deserialization");
        let mut res = Self::default();
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];
        let mut buf = [0u8];

        res.nodes.iter_mut().for_each(|(left, right)| {
            left.iter_mut().for_each(|poly| {
                poly.coeffs.iter_mut().for_each(|coeff| {
                    reader.read_exact(&mut buf4).unwrap();
                    *coeff = i32::from_le_bytes(buf4);
                })
            });
            right
                .iter_mut()
                .take(HVC_WIDTH - skip_last_poly as usize)
                .for_each(|poly| {
                    poly.coeffs.iter_mut().for_each(|coeff| {
                        reader.read_exact(&mut buf4).unwrap();
                        *coeff = i32::from_le_bytes(buf4);
                    })
                });
        });
        reader.read_exact(&mut buf8).unwrap();
        res.index = usize::from_le_bytes(buf8);
        reader.read_exact(&mut buf).unwrap();
        assert!(buf[0] == 0 || buf[0] == 1);
        res.is_randomized = buf[0] != 0;
        end_timer!(timer);
        res
    }
}

#[cfg(test)]
mod test {

    use std::io::Cursor;

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
            {
                // serialization
                let mut bytes = vec![];
                path.serialize(&mut bytes, true, false);
                let buf = Cursor::new(bytes);
                let path_reconstructed = RandomizedPath::deserialize(buf, false);
                assert_eq!(path, path_reconstructed);

                // complete
                let mut path2 = path.clone();
                path2
                    .nodes
                    .iter_mut()
                    .for_each(|(_left, right)| right[HVC_WIDTH - 1] = HVCPoly::default());

                path2
                    .nodes
                    .iter()
                    .enumerate()
                    .for_each(|(i, (_left, right))| println!("{}-th: {}", i, right[HVC_WIDTH - 1]));

                path2.complete(&roots, &hasher);

                path2
                    .nodes
                    .iter()
                    .enumerate()
                    .for_each(|(i, (_left, right))| println!("{}-th: {}", i, right[HVC_WIDTH - 1]));

                // assert_eq!(path, path2);
            }
            assert!(path.verify(&roots, &hasher));
        }
    }
}
