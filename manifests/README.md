# Processing manifests

This directory contains compact, machine-readable manifests used to connect
sample metadata with pipeline inputs.

## `processing_units.tsv`

Unified processing-unit table for Angkor, Nanzuo, and MCP. Each row represents
one independently traceable input unit and records its project, sample/library
identifiers, assay, processing bucket, source path, best-hit readiness, and
read-name tagging format.

Snapshot provenance:

- Source: `/Users/inmpain/Desktop/rice/manifests/processing_units.tsv`
- Snapshot date: 2026-09-10
- Data rows: 1,132
- SHA-256: `ba11e9602efaf6ef6b45ef46661bad6c34731128312e33360dc72d34d36bed4c`

Important caveat: MCP `proxy...` values are temporary processing keys, not
confirmed final biological-sample identifiers. Nanzuo `function` units already
use IRGSP-coordinate BAMs and therefore bypass competitive best-hit mapping.

## `rice_environmental_sample_master_v3.csv`

One row per actual environmental sample, covering Angkor, Nanzuo, MCP, and the
non-Angkor/QC519 records retained in the reviewed master table.

- Snapshot date: 2026-09-11
- Data rows: 806
- Columns: 63
- Unique `sample_id`: 806
- SHA-256: `b614f4eb399c23b9a25b339cd2bfd3ae50d9c52efa07196be5e54e216e0e012e`

Coordinate interpretation:

- `reported_in_cam_metadata` (376): Angkor coordinates copied by exact
  `core_id` from the `field_sample_id@field_sample`, `latitude`, and `longitude`
  columns in `sample_meta_data_20250922.tsv`. The five core assignments were
  independently reproduced in `angkor_final_metadata.tsv`.
- `reported_in_sample_metadata` (33): copied unchanged from the Nanzuo/MCP
  sample metadata.
- `needs_georeferencing` (36): a site name is present but a reliable sampling
  coordinate is not yet available.
- `missing_site_metadata` (361): no site was available, so coordinates remain
  blank rather than being guessed.

Use `coordinate_status`, `coordinate_source`, and `coordinate_note` whenever
filtering or replacing coordinates. Empty latitude/longitude cells are missing
data, not zero coordinates.
