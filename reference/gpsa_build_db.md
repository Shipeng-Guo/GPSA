# Build GPSA HDF5 Database

Constructs an HDF5 database from a diffTable RDS file containing genes x
signatures logFC matrix. The database stores rankZ-transformed values
for fast similarity search.

## Usage

``` r
gpsa_build_db(
  rds_file,
  h5file,
  meta_file = NULL,
  gene_col = NULL,
  na_fill_rank = 0,
  ties_method = "average",
  block_s = 128,
  chunk_g_Z = 1024,
  chunk_s_Z = 512,
  chunk_g_logFC = 2048,
  compress_level = 1,
  store_logFC = TRUE,
  use_float32 = TRUE,
  verbose = TRUE
)
```

## Arguments

- rds_file:

  Path to RDS file containing diffTable (genes x signatures)

- h5file:

  Output HDF5 file path

- meta_file:

  Optional path to metadata RDS (must contain signature_id, pbgene,
  CellLineName)

- gene_col:

  Column name for gene symbols (default: first column)

- na_fill_rank:

  Value to fill NA before ranking (default: 0)

- ties_method:

  Ties method for ranking (default: "average")

- block_s:

  Block size for writing (default: 128)

- chunk_g_Z:

  Chunk size for genes in Z matrix (default: 1024)

- chunk_s_Z:

  Chunk size for signatures in Z matrix (default: 512)

- chunk_g_logFC:

  Chunk size for genes in logFC matrix (default: 2048)

- compress_level:

  Compression level 0-9 (default: 1, lower = faster read)

- store_logFC:

  Store raw logFC values? (default: TRUE)

- use_float32:

  Use float32 storage? (default: TRUE)

- verbose:

  Print progress? (default: TRUE)

## Value

TRUE invisibly on success

## Examples

``` r
if (FALSE) { # \dontrun{
gpsa_build_db(
  rds_file = "diffTable.rds",
  h5file = "gpsa_db.h5",
  meta_file = "metadata.rds"
)
} # }
```
