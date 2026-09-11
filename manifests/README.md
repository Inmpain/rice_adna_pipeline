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
- SHA-256: `5f5602df8efd8b3f0468b58472014e8da2d0ae3d5565850a925d9f0ad0400022`

Coordinate interpretation:

- `reported_in_sample_metadata` (33): copied unchanged from the Nanzuo/MCP
  sample metadata.
- `approximate_site_coordinate` (376): Angkor locality-level coordinates from
  cited gazetteer/map sources; these are not sediment-core GPS coordinates.
- `needs_georeferencing` (36): a site name is present but a reliable sampling
  coordinate is not yet available.
- `missing_site_metadata` (361): no site was available, so coordinates remain
  blank rather than being guessed.

Use `coordinate_status`, `coordinate_source`, and `coordinate_note` whenever
filtering or replacing coordinates. Empty latitude/longitude cells are missing
data, not zero coordinates.
