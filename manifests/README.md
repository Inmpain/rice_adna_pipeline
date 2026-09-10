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

