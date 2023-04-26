use rand::Rng;
use crate::Polynomial;
use crate::{HOTSNTTPoly, GAMMA};

// HOTS public parameters
#[derive(Debug, Clone, Copy)]
pub struct HOTSParam {
    pub(crate) a: [HOTSNTTPoly; GAMMA],
}

impl HOTSParam {
    pub(crate) fn setup<R: Rng>(rng: &mut R) -> Self {
        let mut a = [HOTSNTTPoly::default(); GAMMA];
        a.iter_mut()
            .for_each(|x| *x = HOTSNTTPoly::from(&crate::HOTSPoly::rand_poly(rng)));

        Self { a }
    }
}
