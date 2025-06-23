#!/bin/bash

# Get the absolute path of the script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Initializing installation with Conda..."

# Verify conda installation
if ! command -v conda &> /dev/null; then
    echo "Conda not found, installing..."
    # Install miniconda if not present
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    source "$HOME/miniconda/bin/activate"
fi

# Update packages and install dependencies via APT
echo "Installing dependencies..."
sudo apt update -y
sudo apt install -y build-essential uuid-dev libgpgme-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup-bin r-base

# Verify env 'poly_pipeline'
if conda env list | grep -q "poly_pipeline"; then
    echo "Environment found, removing..."
    conda deactivate
    conda env remove --name poly_pipeline -y || { echo "Failed during environment removal, exiting."; exit 1; }
fi

# Create env 'poly_pipeline'
echo "Creating Conda Environment"
conda env create -f environment.yml -y || { echo "Failed during environment creation, exiting."; exit 1; }

# Install R packages
echo "Installing R packages..."
conda run -n poly_pipeline R -e "install.packages('ggplot2')"

#!/bin/bash
# Compile C++ program in bin/
echo "Compiling executable..."
mkdir -p "$SCRIPT_DIR/bin"

# Check if main.cpp exists
if [ ! -f "$SCRIPT_DIR/main.cpp" ]; then
    echo "Error: main.cpp not found in $SCRIPT_DIR"
    exit 1
fi

echo "Compiling $SCRIPT_DIR/main.cpp..."
g++ "$SCRIPT_DIR/main.cpp" -o "$SCRIPT_DIR/bin/poly_exe" `pkg-config --cflags --libs gtkmm-3.0` || {
    echo "Compilation failed!"
    exit 1
}

# Generate script poly_app
cat << EOF > "$SCRIPT_DIR/poly_app"
#!/bin/bash
# Activate Conda and initiate the program

SCRIPT_DIR="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )"

# Ensure the 'inputs' directory exists
mkdir -p "$SCRIPT_DIR/INPUT"

# Ensure the 'inputs' directory exists
mkdir -p "$SCRIPT_DIR/RESULTS"

eval "\$(conda shell.bash hook)"
conda activate poly_pipeline
exec "\$SCRIPT_DIR/bin/poly_exe" "\$@"
EOF

chmod +x "$SCRIPT_DIR/poly_app"

echo "✅ All set! Run ./poly_app to run the program!"
