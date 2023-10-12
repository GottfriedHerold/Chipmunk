//! This module implements serialization and the trick to reduce ring element by 1
//!

use std::io::{Read, Write};

use ark_std::{end_timer, start_timer};

use crate::{randomize_path::RandomizedPath, ChipmunkSignature, HotsSig, RandomizedHOTSPK};

impl ChipmunkSignature {
    pub(crate) fn serialize<W: Write>(
        &self,
        mut writer: W,
        is_aggregated: bool,
        skip_last_poly: bool,
    ) {
        let timer = start_timer!(|| "Chipmunk serialization");
        self.path
            .serialize(&mut writer, is_aggregated, skip_last_poly);
        self.hots_pk.serialize(&mut writer, is_aggregated);
        self.hots_sig.serialize(&mut writer, is_aggregated);
        end_timer!(timer);
    }

    pub(crate) fn deserialize<R: Read>(mut reader: R, skip_last_poly: bool) -> Self {
        let timer = start_timer!(|| "Chipmunk deserialization");
        let mut res = Self::default();
        res.path = RandomizedPath::deserialize(&mut reader, skip_last_poly);
        res.hots_pk = RandomizedHOTSPK::deserialize(&mut reader);
        res.hots_sig = HotsSig::deserialize(&mut reader);

        end_timer!(timer);
        res
    }
}
