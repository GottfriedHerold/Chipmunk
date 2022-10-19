#!/bin/bash -e

export RUST_MIN_STACK=100000000
echo "      date     time $(free -g | grep total | sed -E 's/^    (.*)/\1/g')"
echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -g | grep Mem: | sed 's/Mem://g')"
cargo test --release benchmark_chipmunk -- --nocapture &
while sleep 10; do
	echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -g | grep Mem: | sed 's/Mem://g')"
done

RAYON_NUM_THREADS=64 cargo test --release benchmark_chipmunk -- --nocapture 
RAYON_NUM_THREADS=32 cargo test --release benchmark_chipmunk -- --nocapture 
RAYON_NUM_THREADS=24 cargo test --release benchmark_chipmunk -- --nocapture 