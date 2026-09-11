
## 2026-07-29 进展更新

完成9格提取方法x定量比对工具矩阵测试，详见docs/09_extraction_mapping_matrix_final.md
和docs/decisions_log.md。核心结论：提取阶段(阶段①)用BWA是决定性因素，
定量比对阶段用BWA或Bowtie2新参数(-N1)均可接受。为"坚持用BWA"这一决策
提供了最完整的证据支撑。

## 2026-09-10 metadata snapshot

Added `manifests/processing_units.tsv`, the compact unified processing-unit
manifest for Angkor, Nanzuo, and MCP (1,132 units). Its provenance and usage
constraints are documented in `manifests/README.md` and
`docs/PROCESSING_UNITS_METADATA_HANDOFF.md`.

## 2026-09-11 environmental sample master v3

Added `manifests/rice_environmental_sample_master_v3.csv`, a one-row-per-sample
table for 806 environmental samples. It preserves the staged read counts from
the reviewed v2 workbook and adds coordinate status, source, and notes. The
table distinguishes reported coordinates, approximate Angkor locality points,
sites still needing georeferencing, and rows lacking site metadata.
