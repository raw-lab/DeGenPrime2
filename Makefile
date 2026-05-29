# DeGenPrime 2.1 Makefile - cross-platform conveniences
#
# Usage:
#   make install    # full install (Rust + Python venv + launchers + icon)
#   make build      # build the Rust binary
#   make test       # run Rust unit tests
#   make gui        # launch the Python GUI
#   make clean      # remove build artifacts
#   make icon       # regenerate the icon assets

.PHONY: help install build test test-verbose gui cli clean rebuild format lint icon

help:
	@echo "DeGenPrime 2.1 - convenience commands"
	@echo ""
	@echo "  make install     Run the full installer (install.sh / install.bat)"
	@echo "  make build       Build the Rust binary in release mode"
	@echo "  make test        Run all Rust unit tests"
	@echo "  make test-verbose"
	@echo "                   Run tests with full output"
	@echo "  make gui         Launch the GUI (after install)"
	@echo "  make cli         Print CLI help"
	@echo "  make icon        Regenerate the desktop icons"
	@echo "  make lint        cargo clippy"
	@echo "  make format      cargo fmt"
	@echo "  make clean       Remove build artifacts"
	@echo "  make rebuild     Clean + build"

install:
	@if [ "$$(uname -s)" = "Linux" ] || [ "$$(uname -s)" = "Darwin" ]; then \
		bash install.sh; \
	else \
		cmd /c install.bat; \
	fi

build:
	cargo build --release

test:
	cargo test --release

test-verbose:
	cargo test --release -- --nocapture

gui:
	@if [ -x run_gui.sh ]; then ./run_gui.sh; else \
	  echo "Run 'make install' first."; exit 1; fi

cli:
	@if [ -x target/release/degenprime ]; then ./target/release/degenprime --help; else \
	  echo "Run 'make build' first."; exit 1; fi

lint:
	cargo clippy --release --all-targets -- -D warnings

format:
	cargo fmt

clean:
	cargo clean
	rm -rf venv build
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

rebuild: clean build

icon:
	python3 scripts/make_icon.py

.DEFAULT_GOAL := help
