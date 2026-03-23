[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# TELE-Vir Installation

This repository contains the installation pipeline for **TELE-Vir** - a tool for the identification of viral sequences in metagenomic data.

**Purpose:** This is an **installation-only** repository. It handles the installation of:

- Conda environments for bioinformatics tools
- Reference databases (NCBI RefSeq, Kraken2, Centrifuge, etc.)
- Host genome sequences
- Software indexes

For the main TELE-Vir analysis pipeline, see the INSaFLU project.

---

## What This Repository Installs

### Bioinformatics Tools

| Category             | Tools                                                                               |
| -------------------- | ----------------------------------------------------------------------------------- |
| Quality Control      | FastQC, Trimmomatic, RabbitQC                                                       |
| Host Depletion       | BWA, Minimap2, Bowtie2                                                              |
| Classification       | Centrifuge, Kraken2, KrakenUniq, Kaiju, FastViromeExplorer, CLARK, deSAMBA, Voyager |
| Remapping            | Bowtie2, Minimap2                                                                   |
| RNA Quantification   | Kallisto                                                                            |
| Microbiome Profiling | MetaPhlAn                                                                           |

### Databases

| Database        | Description                                                                                                                                   |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| RefSeq Genomes  | NCBI viral/bacterial genome sequences                                                                                                         |
| RefSeq Proteins | Protein sequences                                                                                                                             |
| 16S rRNA        | RefSeq, SILVA, NCBI RDP databases                                                                                                             |
| Protein         | UniRef90, UniRef100, SwissProt, RVDB, Virosaurus                                                                                              |
| Host Genomes    | 15 species: human (hg38), cat, dog, pig, cow, chicken, duck, salmon, rainbow trout, mink, marmot, bat, mosquito (Culex, Aedes), sandfly, carp |
| Kraken2         | viral, standard, bacteria_16gb, RDP 16S, eupathdb48                                                                                           |
| Centrifuge      | viral, bacteria                                                                                                                               |
| MetaPhlAn       | CHOCOPhlAn SGB (default, vJan25)                                                                                                              |
| Kaiju           | viruses_2024, fungi_2024, refseq_2024                                                                                                         |
| Voyager         | viral, bacteria                                                                                                                               |
| Taxonomy        | NCBI taxdump, accession2taxid                                                                                                                 |

---

## Quick Start

### Build the Docker Image

```bash
docker build -t televir .
```

### Run Installation

```bash
# Installation with persistence to /data
docker run -v /data:/opt/televir televir move
```

---

## Directory Structure

The container has two separate locations:

| Path                | Purpose                                  | Persisted |
| ------------------- | ---------------------------------------- | --------- |
| `/opt/televir-repo` | Repository code (read-only)              | No        |
| `/opt/televir`      | Data directory (environments, databases) | Yes       |

### Inside the Container (Image)

```
/opt/televir-repo/        # Repository code (downloaded at build time)
├── install_scripts/      # Installation scripts
├── main.py              # Entry point
├── yaml/                # Conda environment definitions
└── ...

/opt/televir/            # Data directory (can be mounted)
├── environments/        # Conda environments
├── ref_db/              # Software database indexes
├── ref_fasta/           # Reference sequences
├── metadata/            # NCBI metadata
├── utility_local.db     # Registration database
├── software.tsv         # Software registry
└── database.tsv         # Database registry
```

### Why Separate Repository from Data?

- **Repository (`/opt/televir-repo`)**: Contains code/scripts, stays in the image
- **Data (`/opt/televir`)**: Contains environments and databases, persisted to volume

This allows you to mount `/opt/televir` to a host directory without polluting it with repository files.

## Entrypoint Commands

This container provides the following commands for installation management:

### `move` - Install with Persistence

```bash
docker run -v /data:/opt/televir televir move
```

**Purpose:** Install TELE-Vir with data persisted to the mounted volume.

**Use case:** First-time installation where you want databases and environments persisted outside the container.

**What it does:**

1. Creates directories in the mounted volume (`/opt/televir`)
2. Installs conda environments
3. Downloads and builds databases
4. Registers installed items in SQLite database
5. Sets appropriate permissions

**Note:** The repository code stays in `/opt/televir-repo` (inside the image) and is not copied to the volume.

---

### `install` - In-Place Installation

```bash
docker run televir install
```

**Purpose:** Install TELE-Vir with data stored in the container's `/opt/televir`.

**Use case:** Testing or temporary installations.

**Warning:** All data will be lost when the container is removed.

---

### `update` - Rebuild Databases

```bash
docker run -v /data:/opt/televir -e UPDATE=true televir update
```

**Purpose:** Update existing installation by rebuilding all databases.

**Use case:** Periodic database updates to get latest RefSeq versions.

**What it does:**

1. Removes existing database directories
2. Re-downloads all databases from NCBI
3. Re-builds all software indexes
4. Re-registers everything with new version information

---

### `check` - Verify Installation

```bash
docker run -v /data:/opt/televir televir check
```

**Purpose:** Display current installation status.

**Use case:** Verify installation completed correctly or check what's installed.

**Output includes:**

- Repository information
- List of installed databases
- List of installed environments
- Registration database status
- Configuration values
- Software versions
- Database versions

---

### `sources` - Query Source Configuration

```bash
# List all sources
docker run televir sources list all

# List databases
docker run televir sources list databases
docker run televir sources list databases -u    # with URLs

# List host genomes
docker run televir sources list hosts
docker run televir sources list hosts -u       # with URLs

# List software
docker run televir sources list software
docker run televir sources list software -u   # with URLs

# Get specific source
docker run televir sources get database kraken2 viral
docker run televir sources get host homo_sapiens
docker run televir sources get software trimmomatic

# Get URL format only
docker run televir sources get host homo_sapiens -f url

# Validate configuration
docker run televir sources validate

# Query with volume mounted (post-install)
docker run -v /data:/opt/televir televir sources list databases
```

**Purpose:** Query and validate the source configuration (databases, software, host genomes).

**Use case:** Inspect available sources, verify URLs, or validate configuration before installation.

**Available after installation:** Tools are copied to `$INSTALL_HOME/tools/`:

```bash
# From mounted volume
docker run -v /data:/opt/televir --entrypoint python televir \
    /opt/televir/tools/sources_cli.py list databases
```

### `status` - GUI Status Display

```bash
# install x-11
sudo apt-get install x11-xserver-utils

# Allow connections to the X server
xhost +local:docker

# Display GUI status
docker run -v /data:/opt/televir -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$DISPLAY televir status
```

**Purpose:** Display GUI showing all databases and software with installation status.

**Features:**

- Two panels showing Databases and Software tables
- Columns: Name, Category/Type, Available, Installed
- "Check Status" button refreshes installation status from utility_local.db

**Note:** Requires graphical display environment or X11 forwarding.

---

## Using Custom Local Directories

By default, TELE-Vir uses `/opt/televir` inside the container. You can customize this:

### Persist to Custom Directory

```bash
# Mount your desired directory
docker run -v /my/custom/path:/opt/televir televir move
```

### Environment Variables

| Variable           | Default                             | Description                 |
| ------------------ | ----------------------------------- | --------------------------- |
| `INSTALL_HOME`     | `/opt/televir`                      | Main installation directory |
| `ENVDIR`           | `/opt/televir/environments`         | Conda environments location |
| `SOURCE`           | `/opt/conda/etc/profile.d/conda.sh` | Conda activation script     |
| `TAXDUMP`          | `/opt/taxdump.tar.gz`               | NCBI taxonomy dump          |
| `UPDATE`           | `false`                             | Set to `true` for rebuild   |
| `REQUEST_SEQ_FILE` | -                                   | Custom request sequences    |

### Examples

```bash
# Persist to /data on host
docker run -v /data:/opt/televir televir move

# Custom data directory
docker run -v /my/televir:/opt/televir televir move

# Update existing installation
docker run -v /data:/opt/televir -e UPDATE=true televir update

# Check installation status
docker run -v /data:/opt/televir televir check

# With custom taxdump
docker run \
  -v /data:/opt/televir \
  -v /path/to/taxdump.tar.gz:/opt/taxdump.tar.gz \
  televir move

# With request sequences
docker run \
  -v /data:/opt/televir \
  -v /path/to/sequences.fa.gz:/data/request_sequences.fa.gz \
  televir move

# Query source configuration
docker run televir sources list all
docker run televir sources list databases -u
docker run televir sources get database kraken2 viral
docker run televir sources validate

# Query from installed environment
docker run -v /data:/opt/televir televir sources list databases

# Display GUI status
docker run -v /data:/opt/televir televir status
```

---

## Request Sequences

You can provide custom request sequences in three ways, checked in this priority order:

### 1. Environment Variable (highest priority)

```bash
docker run -e REQUEST_SEQ_FILE=/data/request_sequences.fa.gz -v /data:/opt/televir televir move
```

### 2. Build-Time (via `--build-arg`)

```bash
docker build --build-arg REQUEST_SEQ_FILE=/path/to/sequences.fa.gz -t televir .
```

### 3. Runtime Volume Mount

```bash
docker run -v /path/to/sequences.fa.gz:/data/request_sequences.fa.gz -v /data:/opt/televir televir move
```

**Priority order:** Environment variable → `/opt/request_sequences.fa.gz` (built-in) → `/data/request_sequences.fa.gz` (mounted)

---

## Directory Structure

After installation, the following structure is created in `/opt/televir` (the data volume):

```
/opt/televir/
├── ref_db/              # Software database indexes
│   ├── kraken2/
│   ├── centrifuge/
│   └── ...
├── ref_fasta/           # Reference sequences
│   ├── refseq_viral.genome.fna.gz
│   └── ...
├── metadata/            # NCBI metadata
├── environments/        # Conda environments
├── tools/              # CLI tools (sources_cli.py, load_sources.py)
├── utility_local.db     # Registration database
├── software.tsv        # Software registry
├── database.tsv        # Database registry
├── config.py           # Runtime configuration
├── install_scripts_config.py  # Installation scripts configuration
└── televir.env         # Environment config
```

**Note:** The repository code is stored separately in `/opt/televir-repo` inside the image and is not copied to the data volume.

---

## Database Version Tracking

TELE-Vir tracks database versions in the registration system:

- **Database table:** Stores name, path, version, source URL, file modification date
- **Software table:** Stores name, path, associated database version, environment path

Version information is captured:

- **RefSeq:** File modification date
- **Hosts:** GCF assembly version
- **Kraken2:** Version from URL (e.g., `20250402`)

---

## Building with Custom Request Sequences

To build the image with request sequences included at build time:

```bash
docker build --build-arg REQUEST_SEQ_FILE=/path/to/sequences.fa.gz -t televir .
```

The file will be copied to `/opt/request_sequences.fa.gz` in the image.

---

## Troubleshooting

### Check Installation Status

```bash
docker run -v /data:/opt/televir televir check
```

### View Logs

```bash
# Run interactively to see all output
docker run -it -v /data:/opt/televir televir move
```

### Common Issues

**Database download failures:**

- Check internet connectivity
- Verify TAXDUMP is accessible

**Permission errors:**

- Ensure the volume mount is writable
- Try running without the `move` command first

---

## Local Installation (Non-Docker)

For local installation (without Docker), see below.

### Dependencies

1. Install Miniconda:

   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. Install the `mngsbench_install.yml` environment:
   ```bash
   conda env create -f mngsbench_install.yml
   conda activate mngsbench_install
   ```

### Configuration

Edit `config.py` in the root directory:

- `SOURCE` - Path to conda.sh
- `HOME` - Root directory for sequence and software database installation
- `ENVDIR` - Directory where environments will be installed
- `TAXDUMP` - Path to ncbi's taxdump.tar.gz file

### Deployment

The `main.py` script accepts these boolean flags:

- `--envs` - Install environments
- `--seqdl` - Download reference sequence databases
- `--soft` - Generate software databases

```bash
# Install environments and databases
python main.py --envs --seqdl

# Install everything
python main.py
```

Run `python main.py -h` for details.

---

## Related Projects

- [INSaFLU](https://github.com/INSaFLU/INSaFLU) - Main analysis pipeline
- [INSaFLU Docker](https://github.com/INSaFLU/docker) - Docker deployment

---

## License

GPL v2

## Authors

João Dourado Santos, Miguel Pinheiro, Daniel Sobral, Joana Isidro, Miguel Pinto, João Paulo Gomes and Vítor Borges

## Funding

INSaFLU development is being co-funded by the European Commission on behalf of OneHealth EJP TELE-Vir project.
https://onehealthejp.eu/jrp-tele-vir/
