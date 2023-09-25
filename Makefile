run:
	cargo run

clippy:
	cargo clippy

test:
	cargo test -- --nocapture --test-threads=1
