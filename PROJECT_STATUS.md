
## 2026-07-29 进展更新

完成9格提取方法x定量比对工具矩阵测试，详见docs/09_extraction_mapping_matrix_final.md
和docs/decisions_log.md。核心结论：提取阶段(阶段①)用BWA是决定性因素，
定量比对阶段用BWA或Bowtie2新参数(-N1)均可接受。为"坚持用BWA"这一决策
提供了最完整的证据支撑。
