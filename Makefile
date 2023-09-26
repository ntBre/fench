run:
	cargo run

clippy:
	cargo clippy

testflags = -- --nocapture

ifdef ARGS
    testflags += --test-threads=1
endif

test:
	cargo test $(testflags) $(ARGS)
