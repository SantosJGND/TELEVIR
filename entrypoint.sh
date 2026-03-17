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
    # Copy config from install_scripts to root
    if [ -f "$INSTALL_HOME/install_scripts/config.py" ]; then
        cp "$INSTALL_HOME/install_scripts/config.py" "$INSTALL_HOME/config.py"
    fi
    
    # Copy televir.env for reference
    if [ -f "/opt/televir/televir.env" ]; then
        cp /opt/televir/televir.env "$INSTALL_HOME/televir.env"
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
    
    if [ "$move_flag" = "true" ]; then
        # Copy TELE-Vir files to target directory
        cp -r /opt/televir/TELEVIR/* "$INSTALL_HOME/"
        cd "$INSTALL_HOME"
        
        copy_request_sequences
        setup_config
    else
        cd /opt/televir/TELEVIR
        copy_request_sequences
    fi
    
    # Run installation
    /opt/venv/bin/python main.py \
        --docker \
        --envs \
        --setup_conda \
        --seqdl \
        --soft \
        --partial
    
    if [ "$move_flag" = "true" ]; then
        set_permissions
    fi
    
    echo "---> Installation complete!"
    echo "---> Installed to: $INSTALL_HOME"
}

# Function to update installation
run_update() {
    echo "---> Updating TELE-Vir installation..."
    echo "     INSTALL_HOME: $INSTALL_HOME"
    echo "     ENVDIR: $ENVDIR"
    
    UPDATE=true
    
    cd "$INSTALL_HOME"
    
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
    
    cd "$INSTALL_HOME"
    
    echo ""
    echo "=== Installation Status ==="
    echo ""
    
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
    
    cd "$INSTALL_HOME"
    
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
        echo "---> Command: move (install + move to custom directory)"
        echo "     This will install TELE-Vir to INSTALL_HOME and persist to mounted volume"
        run_install "true"
        ;;
    install)
        echo "---> Command: install (in-place installation)"
        echo "     This will install TELE-Vir in the container's /opt/televir"
        run_install "false"
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
        echo "  move     - Install TELE-Vir and persist to INSTALL_HOME (use with volume mount)"
        echo "  install  - Install TELE-Vir in place (container internal storage)"
        echo "  update   - Update existing installation (rebuild all databases)"
        echo "  check    - Verify installation status"
        echo "  register - Re-register databases (after manual changes)"
        echo ""
        echo "Environment variables:"
        echo "  INSTALL_HOME      - Installation directory (default: /opt/televir)"
        echo "  ENVDIR            - Environments directory (default: /opt/televir/environments)"
        echo "  SOURCE            - Conda source script (default: /opt/conda/etc/profile.d/conda.sh)"
        echo "  TAXDUMP           - Taxdump file location (default: /opt/taxdump.tar.gz)"
        echo "  UPDATE            - Update mode (default: false)"
        echo "  REQUEST_SEQ_FILE  - Optional request sequences file"
        echo ""
        echo "Examples:"
        echo "  # Persist to /data on host"
        echo "  docker run -v /data:/data televir move"
        echo ""
        echo "  # Persist to custom directory"
        echo "  docker run -v /my/televir:/opt/televir televir move"
        echo ""
        echo "  # Update existing installation"
        docker run -v /data:/data -e UPDATE=true televir update
        echo ""
        echo "  # Check installation status"
        echo "  docker run -v /data:/data televir check"
        exit 1
        ;;
esac
