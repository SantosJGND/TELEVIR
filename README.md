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

| Category | Tools |
|----------|-------|
| Quality Control | FastQC, Trimmomatic, RabbitQC |
| Host Depletion | BWA, Minimap2, Bowtie2 |
| Classification | Centrifuge, Kraken2, KrakenUniq, Kaiju, FastViromeExplorer |
| Remapping | Bowtie2, Minimap2 |

### Databases

| Database | Description |
|----------|-------------|
| RefSeq Genomes | NCBI viral/bacterial genome sequences |
| RefSeq Proteins | Protein sequences |
| Host Genomes | 15 host species for depletion |
| Kraken2 Indexes | Pre-built Kraken2 databases |
| Centrifuge Indexes | Centrifuge classification indexes |

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

| Path | Purpose | Persisted |
|------|---------|-----------|
| `/opt/televir-repo` | Repository code (read-only) | No |
| `/opt/televir` | Data directory (environments, databases) | Yes |

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

This container provides five commands for installation management:

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

### `register` - Re-register Databases

```bash
docker run -v /data:/opt/televir televir register
```

**Purpose:** Re-register databases without rebuilding.

**Use case:** After manually adding or modifying database files.

---

## Using Custom Local Directories

By default, TELE-Vir uses `/opt/televir` inside the container. You can customize this:

### Persist to Custom Directory

```bash
# Mount your desired directory
docker run -v /my/custom/path:/opt/televir televir move
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `INSTALL_HOME` | `/opt/televir` | Main installation directory |
| `ENVDIR` | `/opt/televir/environments` | Conda environments location |
| `SOURCE` | `/opt/conda/etc/profile.d/conda.sh` | Conda activation script |
| `TAXDUMP` | `/opt/taxdump.tar.gz` | NCBI taxonomy dump |
| `UPDATE` | `false` | Set to `true` for rebuild |
| `REQUEST_SEQ_FILE` | - | Custom request sequences |

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
```

---

## Request Sequences

You can provide custom request sequences in three ways:

### 1. Build-Time (via `--build-arg`)

```bash
docker build --build-arg REQUEST_SEQ_FILE=/path/to/sequences.fa.gz -t televir .
```

### 2. Runtime Volume Mount

```bash
docker run -v /path/to/sequences.fa.gz:/data/request_sequences.fa.gz -v /data:/opt/televir televir move
```

### 3. Environment Variable

```bash
docker run -e REQUEST_SEQ_FILE=/data/request_sequences.fa.gz -v /data:/opt/televir televir move
```

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
├── utility_local.db     # Registration database
├── software.tsv        # Software registry
├── database.tsv        # Database registry
├── config.py           # Runtime configuration
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
