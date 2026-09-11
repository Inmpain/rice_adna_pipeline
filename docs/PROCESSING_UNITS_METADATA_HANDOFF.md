# Angkor, Nanzuo, and MCP processing-unit metadata

Updated: 2026-09-11

The canonical compact snapshot is `manifests/processing_units.tsv`. It contains
1,132 data rows and uses `unit_id` as the unique processing key. Use this table
for pooling, read-name tagging, unit splitting, and deciding which inputs enter
competitive best-hit mapping.

Key fields are `project`, `biological_sample_id`, `robot_sample_id`,
`library_id`, `assay`, `pool_bucket`, `input_stage`, `input_path_1`,
`input_path_2`, `besthit_ready`, `next_step`, and `tag_format`.

Interpretation constraints:

- A unit is a processing entity; it is not necessarily an independent
  biological sample.
- MCP `proxy...` identifiers remain provisional until a formal sample
  crosswalk is available.
- Nanzuo `function` units are direct-IRGSP inputs and do not enter competitive
  best-hit mapping.
- Controls must remain identifiable through `control_type` and must not be
  silently treated as biological samples.
- Server paths are recorded provenance and are not evidence that a file still
  exists; production runs should perform a read-only path preflight.

When replacing the snapshot, update its row count and SHA-256 in
`manifests/README.md` in the same commit.

## One-row-per-sample table

`manifests/rice_environmental_sample_master_v3.csv` is the companion table for
sample-level analyses. It contains 806 unique `sample_id` rows and retains the
Asian/IRGSP extract, best-hit, mapped, and deduplicated read counts.

The two tables have different grains: `processing_units.tsv` is one row per
processing input, while the v3 CSV is one row per actual sample. Do not join
them by row number. Join with the available `sample_id`, `library_id`, and
`robot_id` identifiers, and keep capture duplicate conflicts unresolved unless
a canonical output is documented.

Coordinate follow-up should prioritize rows where `coordinate_status` is
`needs_georeferencing`, then `missing_site_metadata`. Angkor rows marked
`reported_in_cam_metadata` are core-specific values from the original CAM
sample metadata; do not replace them with locality-level web coordinates.
