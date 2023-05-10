use super::*;
use ark_std::{end_timer, start_timer};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;

#[test]
fn test_chipmunk() {
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
    env_logger::init();

    let mut seed = [0u8; 32];
    let mut rng = ChaCha20Rng::from_seed(seed);

    let pp = Chipmunk::setup(&mut rng);
    let index = 1;
    let size = 10;
    // key gen timer
    {
        let key_gen_timer = start_timer!(|| format!("run key gen {} times", size));
        for _i in 0..size {
            rng.fill_bytes(&mut seed);
            let _ = Chipmunk::key_gen(&seed, &pp);
        }
        end_timer!(key_gen_timer);
    }

    // signing timer
    {
        let (pk, sk) = Chipmunk::key_gen(&seed, &pp);
        let sign_timer = start_timer!(|| format!("run signing {} times", size));
        for i in 0..size {
            let message = format!("this is the {}-th message to sign", i);
            let _ = Chipmunk::sign(&sk, index as usize, message.as_ref(), &pp);
        }
        end_timer!(sign_timer);
    }

    // verification timer
    {
        let (pk, sk) = Chipmunk::key_gen(&seed, &pp);
        let message = "this is the message to sign";
        let sig = Chipmunk::sign(&sk, 0, message.as_ref(), &pp);
        let verify_timer = start_timer!(|| format!("run verify {} times", size));
        for i in 0..size {
            assert!(Chipmunk::verify(&pk, message.as_ref(), &sig, &pp));
        }
        end_timer!(verify_timer);
    }

    // aggregation timer
    {
        let log_size = 7;
        let size = 1 << log_size;

        let mut sigs = vec![];
        let mut pks = vec![];
        let message = "this is the message to sign";
        for _i in 0..size {
            rng.fill_bytes(&mut seed);
            let (pk, sk) = Chipmunk::key_gen(&seed, &pp);
            let sig = Chipmunk::sign(&sk, 0, message.as_ref(), &pp);
            pks.push(pk);
            sigs.push(sig);
            drop(sk);
        }

        for _ in 0..10 - log_size {
            sigs = [sigs.clone(), sigs].concat();
            pks = [pks.clone(), pks].concat();
        }
        println!("=========================================");
        println!("length of pk and sig: {} {}", sigs.len(), pks.len());

        println!("start aggregation");
        let agg_sig = Chipmunk::aggregate(&sigs, &pks);

        println!("start batch verification");
        assert!(Chipmunk::batch_verify(
            &pks,
            message.as_ref(),
            &agg_sig,
            &pp
        ));

        for _ in 0..3 {
            sigs = [sigs.clone(), sigs].concat();
            pks = [pks.clone(), pks].concat();
        }
        println!("=========================================");
        println!("length of pk and sig: {} {}", sigs.len(), pks.len());

        println!("start aggregation");
        let agg_sig = Chipmunk::aggregate(&sigs, &pks);

        println!("start batch verification");
        assert!(Chipmunk::batch_verify(
            &pks,
            message.as_ref(),
            &agg_sig,
            &pp
        ));
    }
}
