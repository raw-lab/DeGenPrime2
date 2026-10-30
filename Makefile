# DeGenPrime — convenience targets
.PHONY: help build test smoke gui doctor wasm web install lint format icon clean
help:
	@echo "make build    release binary (target/release/degenprime)"
	@echo "make test     unit tests"
	@echo "make smoke    end-to-end test of the CLI and the desktop GUI's API"
	@echo "make gui      open the desktop app"
	@echo "make doctor   check the installation"
	@echo "make wasm     rebuild web/assets/degenprime.wasm"
	@echo "make web      package the Netlify site and WordPress plugin into dist/"
	@echo "make install  run install.sh"
	@echo "make icon     rebuild every icon from assets/logo-source.png (needs Pillow and llvm)"
	@echo "make lint | format | clean"
build:
	cargo build --release
test:
	cargo test --release
smoke:
	bash scripts/smoke_test.sh
gui:
	cargo run --release -- gui
doctor:
	cargo run --release -- doctor
wasm:
	bash scripts/build_wasm.sh
web:
	bash scripts/build_web.sh
install:
	bash install.sh
lint:
	cargo clippy --release --all-targets -- -D warnings
format:
	cargo fmt
icon:
	python3 scripts/make_icon.py
	bash scripts/make_windows_resource.sh
clean:
	cargo clean
	rm -rf dist
.DEFAULT_GOAL := help
