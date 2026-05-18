.PHONY: help setup build install test clean gui cli docs

# DeGenPrime 2.0 Makefile

help:
	@echo "DeGenPrime 2.0 - Available Commands"
	@echo ""
	@echo "Setup & Installation:"
	@echo "  make setup          Run setup script (install dependencies)"
	@echo "  make build          Build Rust backend in release mode"
	@echo "  make install        Install Python GUI dependencies"
	@echo ""
	@echo "Running:"
	@echo "  make gui            Launch GUI application"
	@echo "  make cli            Launch CLI (requires input file)"
	@echo ""
	@echo "Development:"
	@echo "  make test           Run all tests"
	@echo "  make test-verbose   Run tests with output"
	@echo "  make docs           Generate and open documentation"
	@echo "  make lint           Check code quality"
	@echo ""
	@echo "Maintenance:"
	@echo "  make clean          Remove build artifacts"
	@echo "  make rebuild        Clean and rebuild"
	@echo "  make format         Format Rust code"
	@echo ""
	@echo "Examples:"
	@echo "  make setup          # Setup everything"
	@echo "  make gui            # Launch GUI"
	@echo ""

setup:
	@echo "Running setup script..."
	@bash setup.sh || cmd /c setup.bat

build:
	@echo "Building Rust backend..."
	cargo build --release
	@echo "✓ Build complete!"

install:
	@echo "Installing Python dependencies..."
	pip3 install -r gui/requirements.txt
	@echo "✓ Installation complete!"

test:
	@echo "Running tests..."
	cargo test --release
	@echo "✓ All tests passed!"

test-verbose:
	@echo "Running tests with output..."
	cargo test --release -- --nocapture --test-threads=1

gui:
	python3 gui/degenprime_gui.py

cli:
	./target/release/degenprime

docs:
	@echo "Generating documentation..."
	cargo doc --no-deps --open

lint:
	@echo "Checking code quality..."
	cargo clippy -- -D warnings

clean:
	@echo "Cleaning build artifacts..."
	cargo clean
	rm -rf target/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete
	@echo "✓ Clean complete!"

rebuild: clean build
	@echo "✓ Rebuild complete!"

format:
	@echo "Formatting code..."
	cargo fmt
	@echo "✓ Format complete!"

check:
	@echo "Checking project..."
	cargo check
	@echo "✓ Check complete!"

bench:
	@echo "Running benchmarks..."
	cargo bench

.DEFAULT_GOAL := help
