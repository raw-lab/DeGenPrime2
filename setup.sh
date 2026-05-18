#!/bin/bash

################################################################################
# DeGenPrime 2.0 - Setup Script
# Automatic installation and configuration for macOS/Linux
################################################################################

set -e

echo "╔════════════════════════════════════════════════════════════╗"
echo "║         DeGenPrime 2.0 - Setup Script                      ║"
echo "║   Fast, Accurate PCR Primer Design Tool                    ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if running on macOS or Linux
OS_TYPE=$(uname -s)

echo -e "${BLUE}[1/7]${NC} Detecting system type..."
if [[ "$OS_TYPE" == "Darwin" ]]; then
    echo -e "${GREEN}✓ macOS detected${NC}"
    IS_MAC=true
elif [[ "$OS_TYPE" == "Linux" ]]; then
    echo -e "${GREEN}✓ Linux detected${NC}"
    IS_MAC=false
else
    echo -e "${RED}✗ Unsupported OS: $OS_TYPE${NC}"
    echo "Supported: macOS, Linux"
    exit 1
fi

echo ""
echo -e "${BLUE}[2/7]${NC} Checking Rust installation..."
if ! command -v rustc &> /dev/null; then
    echo -e "${YELLOW}⚠ Rust not found. Installing...${NC}"
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source $HOME/.cargo/env
    echo -e "${GREEN}✓ Rust installed${NC}"
else
    RUST_VERSION=$(rustc --version)
    echo -e "${GREEN}✓ $RUST_VERSION found${NC}"
fi

echo ""
echo -e "${BLUE}[3/7]${NC} Checking Python installation..."
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}✗ Python 3 not found${NC}"
    if [[ "$IS_MAC" == true ]]; then
        echo "Install with: brew install python3"
    else
        echo "Install with: apt-get install python3"
    fi
    exit 1
else
    PYTHON_VERSION=$(python3 --version)
    echo -e "${GREEN}✓ $PYTHON_VERSION found${NC}"
fi

echo ""
echo -e "${BLUE}[4/7]${NC} Checking for MAFFT (optional)..."
if command -v mafft &> /dev/null; then
    MAFFT_VERSION=$(mafft --version 2>&1 | head -n 1)
    echo -e "${GREEN}✓ MAFFT found: $MAFFT_VERSION${NC}"
else
    echo -e "${YELLOW}⚠ MAFFT not found (optional, built-in fallback available)${NC}"
    read -p "Install MAFFT now? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        if [[ "$IS_MAC" == true ]]; then
            echo "Installing MAFFT via Homebrew..."
            brew install mafft
            echo -e "${GREEN}✓ MAFFT installed${NC}"
        else
            echo "Installing MAFFT via apt..."
            sudo apt-get update
            sudo apt-get install -y mafft
            echo -e "${GREEN}✓ MAFFT installed${NC}"
        fi
    else
        echo -e "${YELLOW}⚠ Skipping MAFFT installation${NC}"
    fi
fi

echo ""
echo -e "${BLUE}[5/7]${NC} Building Rust backend..."
echo "This may take 2-3 minutes on first build..."
cargo build --release
echo -e "${GREEN}✓ Rust backend built${NC}"

echo ""
echo -e "${BLUE}[6/7]${NC} Installing Python GUI dependencies..."
pip3 install --upgrade pip
pip3 install -r gui/requirements.txt
echo -e "${GREEN}✓ Python dependencies installed${NC}"

echo ""
echo -e "${BLUE}[7/7]${NC} Creating launch scripts..."

# Create convenient launch scripts
cat > run_gui.sh << 'EOF'
#!/bin/bash
python3 gui/degenprime_gui.py
EOF
chmod +x run_gui.sh

cat > run_cli.sh << 'EOF'
#!/bin/bash
./target/release/degenprime "$@"
EOF
chmod +x run_cli.sh

echo -e "${GREEN}✓ Launch scripts created${NC}"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              Setup Complete! 🎉                            ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo -e "${GREEN}Quick start:${NC}"
echo ""
echo "  GUI (recommended for most users):"
echo "    ${BLUE}./run_gui.sh${NC}"
echo "    or"
echo "    ${BLUE}python3 gui/degenprime_gui.py${NC}"
echo ""
echo "  Command Line (for power users):"
echo "    ${BLUE}./run_cli.sh --input_file sequences.fasta${NC}"
echo "    or"
echo "    ${BLUE}./target/release/degenprime --input_file sequences.fasta${NC}"
echo ""
echo "  Get help:"
echo "    ${BLUE}./target/release/degenprime --help${NC}"
echo ""
echo -e "${GREEN}Next steps:${NC}"
echo "  1. Read INSTALL.md for detailed documentation"
echo "  2. Check EXAMPLES.md for usage examples"
echo "  3. Review COMPLETE_SUMMARY.md for feature overview"
echo ""
