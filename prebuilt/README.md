# Pre-built Indices

This directory is for user-provided pre-built centrifuge and kraken2 indices.

## Structure

```
prebuilt/
├── centrifuge/
│   ├── viral/       # Place Centrifuge viral index files here
│   └── bacteria/    # Place Centrifuge bacterial index files here
└── kraken2/
    ├── viral/       # Place Kraken2 viral index files here
    ├── bacteria/    # Place Kraken2 bacterial index files here
    └── eupathdb46/  # Place Kraken2 eupathdb46 index files here
```

## Usage

1. Place your pre-built index files in the appropriate subdirectory
2. Optionally modify the path in `sources.yaml` if you use a custom location
3. During installation, the script will detect and use these indices instead of downloading

## Expected Files

- **Centrifuge**: Look for `.1.cf`, `.2.cf`, etc. index files
- **Kraken2**: Look for `taxo.k2d` file in the database directory