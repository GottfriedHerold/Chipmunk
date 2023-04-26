use super::*;
use ark_std::{start_timer, end_timer};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;

#[test]
fn test_chipmunk() {
    env_logger::init();

    let message = "this is the message to sign";
    let mut seed = [0u8; 32];
    let mut rng = ChaCha20Rng::from_seed(seed);

    let pp = Chipmunk::setup(&mut rng);

    for _ in 0..10 {
        rng.fill_bytes(&mut seed);
        let (pk, sk) = Chipmunk::key_gen(&seed, &pp);
        for _ in 0..10 {
            let index = rng.next_u32() % (1 << HEIGHT - 1);
            let sig = Chipmunk::sign(&sk, index as usize, message.as_ref(), &pp);
            assert!(Chipmunk::verify(&pk, message.as_ref(), &sig, &pp))
        }
    }

    for index in 0..10 {
        let mut sigs = Vec::new();
        let mut pks = Vec::new();
        for _ in 0..10 {
            rng.fill_bytes(&mut seed);
            let (pk, sk) = Chipmunk::key_gen(&seed, &pp);

            let sig = Chipmunk::sign(&sk, index as usize, message.as_ref(), &pp);
            assert!(Chipmunk::verify(&pk, message.as_ref(), &sig, &pp));
            pks.push(pk);
            sigs.push(sig);
        }

        let agg_sig = Chipmunk::aggregate(&sigs, &pks);

        assert!(Chipmunk::batch_verify(
            &pks,
            message.as_ref(),
            &agg_sig,
            &pp
        ))
    }
}

#[test]
fn benchmark_chipmunk() {
    let message = "this is the message to sign";
    let mut seed = [0u8; 32];
    let mut rng = ChaCha20Rng::from_seed(seed);

    let pp = Chipmunk::setup(&mut rng);
    let index = 1;
    let mut sigs = Vec::new();
    let mut pks = Vec::new();
    let mut sks = Vec::new();

    let size = 32;

    let key_gen_timer = start_timer!(|| format!("run key gen {} times", size));
    for _i in 0..size {
        rng.fill_bytes(&mut seed);
        let (pk, sk) = Chipmunk::key_gen(&seed, &pp);
        pks.push(pk);
        sks.push(sk);
    }
    end_timer!(key_gen_timer);

    let sign_timer = start_timer!(|| format!("run signing {} times", size));
    for i in 0..size {
        let sig = Chipmunk::sign(&sks[i], index as usize, message.as_ref(), &pp);
        sigs.push(sig);
    }
    end_timer!(sign_timer);

    let verify_timer = start_timer!(|| format!("run verify {} times", size));
    for i in 0..size {
        assert!(Chipmunk::verify(&pks[i], message.as_ref(), &sigs[i], &pp));
    }
    end_timer!(verify_timer);

    for _ in 0..5 {
        sigs = [sigs.clone(), sigs].concat();
        pks = [pks.clone(), pks].concat();
    }

    println!("length of pk and sig: {} {}", sigs.len(), pks.len());

    println!("start aggregation");
    let agg_sig = Chipmunk::aggregate(&sigs, &pks);

    println!("start batch verification");
    assert!(Chipmunk::batch_verify(
        &pks,
        message.as_ref(),
        &agg_sig,
        &pp
    ))
}
