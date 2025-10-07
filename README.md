# phylosql

Tools for working with microbiome count data stored in an MTGN-style MariaDB
schema. The package focuses on turning the database tables that describe sample
metadata, taxonomy assignments, and sequence variant counts into ready-to-use
[`phyloseq`](https://joey711.github.io/phyloseq/) objects.

## Features

- Opinionated helpers for establishing and caching database connections
- Wrappers that fetch sequence variant counts together with curated taxonomy
  metadata
- Utilities for aggregating ("glomming") abundances to higher taxonomy ranks
  before constructing sparse `phyloseq` objects

## Getting started

1. Install the package and its dependencies (for example with
   `remotes::install_github("your-org/phylosql")`).
2. Create a credentials CSV that contains `host`, `dbname`, `port`, and `user`
   columns and keep the associated password or token handy.
3. Establish a connection and cache it for reuse:

   ```r
   library(phylosql)

   con <- get_mtgn_connection(path = "~/mtgn-creds.csv", key = Sys.getenv("MTGN_KEY"))
   ```

The cached connection can be retrieved implicitly by downstream helpers, so you
rarely need to pass `con` explicitly once it is set.

## Building phyloseq objects

The package exposes a family of functions that collect counts, taxonomy, and
sample metadata and return `phyloseq` objects built on sparse matrices.

- `sql_phyloseq_by_sample()` fetches the full set of sequence variants (SVs) for
  the requested `taxa_group` and optional `samples` filter.
- `sql_phyloseq_by_tax_glom()` aggregates SV counts at a higher taxonomy level
  (for example `"Genus"`) before constructing the object.
- `sql_phyloseq_by_sample_and_tax_glom()` first restricts the dataset to
  selected `samples` and then gloms to the requested taxonomy level.

Each helper validates the requested tables and columns, ensures that taxonomy and
sample metadata exist for the returned SVs, and surfaces clear error messages
when inputs are inconsistent.

```r
phy <- sql_phyloseq_by_tax_glom(
  taxalevel = "Genus",
  taxa_group = "bacteria"
)
```

## Development

- Unit tests can be executed locally with `testthat::test_local()`.
- Database dependent tests require access to an MTGN instance populated with the
  expected tables.

Bug reports and contributions are welcome via pull requests.
