#!/bin/bash

set -e

# Load environment variables from image
if [ -f "/opt/televir/televir.env" ]; then
    source /opt/televir/televir.env
fi

# Allow override from environment
INSTALL_HOME=${INSTALL_HOME:-/opt/televir}
ENVDIR=${ENVDIR:-/opt/televir/environments}
SOURCE=${SOURCE:-/opt/conda/etc/profile.d/conda.sh}
TAXDUMP=${TAXDUMP:-/opt/taxdump.tar.gz}
UPDATE=${UPDATE:-false}

# Request sequences file - check multiple locations (priority: env > build arg > volume)
REQUEST_SEQ_FILE=""
if [ -n "${REQUEST_SEQ_FILE}" ]; then
    : # Use provided value
    elif [ -f "/opt/request_sequences.fa.gz" ]; then
    REQUEST_SEQ_FILE="/opt/request_sequences.fa.gz"
    elif [ -f "/data/request_sequences.fa.gz" ]; then
    REQUEST_SEQ_FILE="/data/request_sequences.fa.gz"
fi

# Function to ensure directories exist
ensure_directories() {
    mkdir -p "$INSTALL_HOME"
    mkdir -p "$ENVDIR"
    mkdir -p "$INSTALL_HOME/ref_db"
    mkdir -p "$INSTALL_HOME/ref_fasta"
    mkdir -p "$INSTALL_HOME/metadata"
}

# Function to copy request sequences if available
copy_request_sequences() {
    if [ -n "$REQUEST_SEQ_FILE" ] && [ -f "$REQUEST_SEQ_FILE" ]; then
        echo "---> Copying request sequences: $REQUEST_SEQ_FILE"
        mkdir -p "$INSTALL_HOME/ref_fasta"
        cp "$REQUEST_SEQ_FILE" "$INSTALL_HOME/ref_fasta/request_references.fa.gz"
    fi
}

# Function to setup config in target directory
setup_config() {
    # Copy root config.py (with paths) to data directory
    if [ -f "/opt/televir-repo/config.py" ]; then
        cp /opt/televir-repo/config.py "$INSTALL_HOME/config.py"
    fi
    
    # Copy config to repo directory for main.py import
    if [ -f "/opt/televir-repo/config.py" ]; then
        cp /opt/televir-repo/config.py /opt/televir-repo/config.py
    fi
    
    # Copy TelevirLayout config to data directory
    if [ -f "/opt/televir-repo/install_scripts/config.py" ]; then
        cp /opt/televir-repo/install_scripts/config.py "$INSTALL_HOME/install_scripts_config.py"
    fi
    
    # Copy televir.env for reference
    if [ -f "/opt/televir-repo/televir.env" ]; then
        cp /opt/televir-repo/televir.env "$INSTALL_HOME/televir.env"
    fi
}

# Function to set permissions
set_permissions() {
    chmod -R 0777 "$INSTALL_HOME"
}

# Function to run installation
run_install() {
    local move_flag="$1"
    
    echo "---> Starting TELE-Vir installation..."
    echo "     INSTALL_HOME: $INSTALL_HOME"
    echo "     ENVDIR: $ENVDIR"
    echo "     UPDATE: $UPDATE"
    
    # Ensure directories exist
    ensure_directories
    
    setup_config
    copy_request_sequences
    
    # Run installation from repo, but install data to /opt/televir
    cd /opt/televir-repo
    /opt/venv/bin/python main.py \
    --docker \
    --envs \
    --setup_conda \
    --seqdl \
    --soft \
    --partial
    
    set_permissions
    
    echo "---> Installation complete!"
    echo "---> Installed to: $INSTALL_HOME"
}

# Function to update installation
run_update() {
    echo "---> Updating TELE-Vir installation..."
    echo "     INSTALL_HOME: $INSTALL_HOME"
    echo "     ENVDIR: $ENVDIR"
    
    UPDATE=true
    cd /opt/televir-repo
    
    copy_request_sequences
    
    /opt/venv/bin/python main.py \
    --docker \
    --envs \
    --setup_conda \
    --seqdl \
    --soft \
    --partial
    
    echo "---> Update complete!"
}

# Function to check installation status
run_check() {
    echo "---> Checking TELE-Vir installation..."
    
    echo ""
    echo "=== TELE-Vir Repository ==="
    echo ""
    if [ -d "/opt/televir-repo" ]; then
        echo "Repository: /opt/televir-repo"
        echo "Contents:"
        ls -la /opt/televir-repo/
    else
        echo "No repository found at /opt/televir-repo"
    fi
    
    echo ""
    echo "=== Installation Status ==="
    
    if [ -d "$INSTALL_HOME/ref_db" ]; then
        echo "Databases installed:"
        ls -la "$INSTALL_HOME/ref_db/"
    else
        echo "No databases found in $INSTALL_HOME/ref_db"
    fi
    
    echo ""
    
    if [ -d "$ENVDIR" ]; then
        echo "Environments installed:"
        ls -la "$ENVDIR/"
    else
        echo "No environments found in $ENVDIR"
    fi
    
    echo ""
    
    if [ -f "$INSTALL_HOME/utility_local.db" ]; then
        echo "Registration database exists"
    else
        echo "No registration database found"
    fi
    
    echo ""
    echo "=== Configuration ==="
    echo ""
    echo "INSTALL_HOME: $INSTALL_HOME"
    echo "ENVDIR: $ENVDIR"
    echo "SOURCE: $SOURCE"
    echo "TAXDUMP: $TAXDUMP"
    
    if [ -n "$REQUEST_SEQ_FILE" ]; then
        echo "REQUEST_SEQ_FILE: $REQUEST_SEQ_FILE"
    fi
    
    echo ""
    echo "=== Software Versions ==="
    echo ""
    
    if [ -f "$INSTALL_HOME/software.tsv" ]; then
        cat "$INSTALL_HOME/software.tsv"
    else
        echo "No software registry found"
    fi
    
    echo ""
    
    if [ -f "$INSTALL_HOME/database.tsv" ]; then
        cat "$INSTALL_HOME/database.tsv"
    else
        echo "No database registry found"
    fi
    
    echo "---> Check complete!"
}

# Function to re-register databases
run_register() {
    echo "---> Re-registering TELE-Vir databases..."
    
    cd /opt/televir-repo
    
    /opt/venv/bin/python main.py \
    --docker \
    --seqdl \
    --soft \
    --partial
    
    echo "---> Registration complete!"
}

# Main command handler
case "$1" in
    move)
        echo "---> Command: move (install to volume mount)"
        echo "     Repository: /opt/televir-repo (in image)"
        echo "     Data: $INSTALL_HOME (mounted volume)"
        run_install
    ;;
    install)
        echo "---> Command: install (in-place installation)"
        echo "     Repository: /opt/televir-repo (in image)"
        echo "     Data: $INSTALL_HOME (internal)"
        run_install
    ;;
    update)
        echo "---> Command: update (rebuild databases)"
        echo "     This will rebuild all databases in the existing installation"
        run_update
    ;;
    check)
        echo "---> Command: check (verify installation)"
        echo "     This will display the current installation status"
        run_check
    ;;
    register)
        echo "---> Command: register (re-register databases)"
        echo "     This will re-register databases (useful after manual DB updates)"
        run_register
    ;;
    *)
        echo "Usage: $0 {move|install|update|check|register} [options]"
        echo ""
        echo "Commands:"
        echo "  move     - Install TELE-Vir (repo in image, data to volume)"
        echo "  install  - Install TELE-Vir (repo in image, data internal)"
        echo "  update   - Update existing installation (rebuild all databases)"
        echo "  check    - Verify installation status"
        echo "  register - Re-register databases (after manual changes)"
        echo ""
        echo "Directory Structure:"
        echo "  /opt/televir-repo  - Repository code (in image, not persisted)"
        echo "  /opt/televir       - Data directory (environments, databases)"
        echo ""
        echo "Environment variables:"
        echo "  INSTALL_HOME      - Data directory (default: /opt/televir)"
        echo "  ENVDIR            - Environments directory (default: /opt/televir/environments)"
        echo "  SOURCE            - Conda source script (default: /opt/conda/etc/profile.d/conda.sh)"
        echo "  TAXDUMP           - Taxdump file location (default: /opt/taxdump.tar.gz)"
        echo "  UPDATE            - Update mode (default: false)"
        echo "  REQUEST_SEQ_FILE  - Optional request sequences file"
        echo ""
        echo "Examples:"
        echo "  # Persist data to /data on host"
        echo "  docker run -v /data:/opt/televir televir move"
        echo ""
        echo "  # Check installation status"
        echo "  docker run -v /data:/opt/televir televir check"
        exit 1
    ;;
esac
