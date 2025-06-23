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
# Install the latest version of BLAST+
echo "Checking for BLAST+ installation..."

# Base URL for BLAST+ downloads
BASE_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"

# Function to extract version number
get_version_number() {
    # Extract version number like 2.16.0+ and remove the '+' character
    echo "$1" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1
}

# Create a temporary directory for downloading the listing
TMP_DIR=$(mktemp -d)
cd "$TMP_DIR" || exit 1

# First, get the latest version number from the NCBI site
echo "Checking latest available BLAST+ version..."
wget -q "$BASE_URL" -O latest_blast_listing.html
if [ ! -f "latest_blast_listing.html" ]; then
    echo "Error: Failed to retrieve the BLAST+ directory listing. Exiting."
    rm -rf "$TMP_DIR"
    exit 1
fi

# Extract the latest version number from the file listing
LATEST_VERSION=$(grep -o "ncbi-blast-[0-9.]\+" latest_blast_listing.html | head -n 1 | sed 's/ncbi-blast-//')
LATEST_VERSION=$(get_version_number "$LATEST_VERSION")
echo "Latest available version: $LATEST_VERSION"

# Variable to control whether we need to install BLAST+
INSTALL_BLAST=true

# Check if BLAST+ is already installed
if command -v blastn &> /dev/null && command -v blastp &> /dev/null; then
    echo "BLAST+ binaries found in PATH."
    # Try to get the version, but handle errors
    CURRENT_VERSION=$(blastn -version 2>&1 | head -n 1)
    
    # Check if there was an error in the version output
    if [[ "$CURRENT_VERSION" == *"error"* ]] || [[ "$CURRENT_VERSION" == *"not found"* ]] || [[ -z "$CURRENT_VERSION" ]]; then
        echo "Found BLAST+ installation appears to be broken or has dependency issues."
        echo "Details of the error:"
        blastn -version 2>&1 | head -n 3
        echo "Will reinstall BLAST+ to fix the issues..."
        
        # Find the location of the broken installation
        BLAST_PATH=$(which blastn)
        echo "Current installation path: $BLAST_PATH"
        
        if [[ -n "$BLAST_PATH" ]]; then
            BLAST_INSTALL_DIR=$(dirname "$(dirname "$BLAST_PATH")")
            echo "Will backup and replace installation at: $BLAST_INSTALL_DIR"
            
            # Backup the old installation if it exists
            if [[ -d "$BLAST_INSTALL_DIR" ]]; then
                echo "Backing up old installation..."
                sudo mv "$BLAST_INSTALL_DIR" "${BLAST_INSTALL_DIR}.bak.$(date +%Y%m%d%H%M%S)"
            fi
        fi
    else
        echo "Current version: $CURRENT_VERSION"
        
        # Extract the version number from the output
        INSTALLED_VERSION=$(get_version_number "$CURRENT_VERSION")
        echo "Installed version number: $INSTALLED_VERSION"
        
        # Compare versions
        if [ "$INSTALLED_VERSION" = "$LATEST_VERSION" ]; then
            echo "You already have the latest version of BLAST+ installed."
            INSTALL_BLAST=false
            
            # Check if BLAST+ is in PATH
            if ! grep -q "/usr/local/blast/bin" ~/.bashrc; then
                echo "Adding BLAST+ to your PATH for future sessions..."
                echo 'export PATH="/usr/local/blast/bin:$PATH"' >> ~/.bashrc
                echo "BLAST+ added to PATH in ~/.bashrc"
                echo "Applying PATH changes to current session..."
                source ~/.bashrc
                echo "PATH updated in current session!"
            fi
        else
            echo "A newer version ($LATEST_VERSION) is available. Proceeding with installation..."
        fi
    fi
else
    echo "BLAST+ not found. Proceeding with installation..."
fi

# Only install BLAST+ if needed
if [ "$INSTALL_BLAST" = true ]; then
    # Determine system architecture
    ARCH=$(uname -m)
    OS=$(uname -s)

    # Set the appropriate file pattern based on architecture and OS
    if [ "$ARCH" = "x86_64" ]; then
        if [ "$OS" = "Linux" ]; then
            FILE_PATTERN="ncbi-blast-.*-x64-linux\.tar\.gz"
        elif [ "$OS" = "Darwin" ]; then
            FILE_PATTERN="ncbi-blast-.*-x64-macosx\.tar\.gz"
        fi
    elif [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then
        if [ "$OS" = "Linux" ]; then
            FILE_PATTERN="ncbi-blast-.*-aarch64-linux\.tar\.gz"
        elif [ "$OS" = "Darwin" ]; then
            FILE_PATTERN="ncbi-blast-.*-aarch64-macosx\.tar\.gz"
        fi
    fi

    if [ -z "$FILE_PATTERN" ]; then
        echo "Unsupported architecture/OS combination: $ARCH/$OS"
        exit 1
    fi

    # Create a temporary directory
    TMP_DIR=$(mktemp -d)
    cd "$TMP_DIR" || exit 1

    # Fetch the latest BLAST+ version
    echo "Finding latest BLAST+ version for $ARCH/$OS..."
    wget -q "$BASE_URL" -O latest_blast_listing.html
    if [ ! -f "latest_blast_listing.html" ]; then
        echo "Error: Failed to retrieve the BLAST+ directory listing. Exiting."
        rm -rf "$TMP_DIR"
        exit 1
    fi

    # Use a more robust approach to find the correct file
    BLAST_FILE=$(grep -o "href=\".*$FILE_PATTERN\"" latest_blast_listing.html | grep -v "\.md5" | head -n 1 | sed 's/href="//;s/"$//')

    if [ -z "$BLAST_FILE" ]; then
        echo "Could not find the download file for pattern: $FILE_PATTERN"
        echo "Available files:"
        grep -o "href=\".*\.tar\.gz\"" latest_blast_listing.html | grep -v "\.md5" | sed 's/href="//;s/"$//'
        rm -rf "$TMP_DIR"
        exit 1
    fi

    echo "Downloading BLAST+ file: $BLAST_FILE"

    # Download and install
    wget "$BASE_URL$BLAST_FILE" -O blast.tar.gz
    if [ ! -f "blast.tar.gz" ]; then
        echo "Failed to download the BLAST+ package. Exiting."
        rm -rf "$TMP_DIR"
        exit 1
    fi

    # Extract to /usr/local
    echo "Installing BLAST+ to /usr/local..."
    sudo mkdir -p /usr/local/blast
    sudo tar -xzf blast.tar.gz -C /usr/local/blast --strip-components=1

    # Cleanup
    cd - || exit 1
    rm -rf "$TMP_DIR"

    # Add to PATH
    echo "Adding BLAST+ to your PATH..."
    if ! grep -q "/usr/local/blast/bin" ~/.bashrc; then
        echo 'export PATH="/usr/local/blast/bin:$PATH"' >> ~/.bashrc
        echo "BLAST+ added to PATH in ~/.bashrc"
    else
        echo "BLAST+ is already in your PATH"
    fi

    echo "BLAST+ installed successfully!"
    echo "Applying PATH changes to current session..."
    source ~/.bashrc
    echo "PATH updated in current session!"
else
    echo "Skipping BLAST+ installation as it's already up-to-date."
fi

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
mkdir -p "$SCRIPT_DIR/inputs"

# Ensure the 'inputs' directory exists
mkdir -p "$SCRIPT_DIR/outputs"

eval "\$(conda shell.bash hook)"
conda activate poly_pipeline
exec "\$SCRIPT_DIR/bin/poly_exe" "\$@"
EOF

chmod +x "$SCRIPT_DIR/poly_app"

echo "✅ All set! Run ./poly_app to run the program!"
