# 更新日志

本文件格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。

## [Unreleased]

## [0.1.0] - 2026-10-23

### 新增

- `makedb` / `blastp` / `blastx` / `dbinfo` 四个子命令，自有 `.dmnd` 风格二进制库（magic `DRDMND01`，v1）。
- 双索引种子检索：DIAMOND 2.2.8 的 spaced seed shape（2×10、16×8、14×7）、`murphy10` / `diamond9`
  约简，以及 `seed_cut * ln2 * weight` 的复杂度过滤。
- ungapped X-drop、对角线分段、重叠 band 合并，以及带状仿射 Smith-Waterman + traceback。
- ALP 有限尺寸校正的 E-value 与 DIAMOND 兼容的 bitscore（`src/alp.rs`）。
- `--outfmt 6` 自定义字段、culling、`--max-target-seqs` / `--max-hsps`。
- 灵敏度预设 `--fast` … `--ultra-sensitive`，`index_chunks` 跟随预设取值。
- `blastx` 六框翻译与核酸坐标映射。
- 基准脚本 `bench/compare.ps1`（对 DIAMOND）与 `bench/compare3.ps1`（再加 NCBI BLAST+）。

### 已知限制

- 未移植 legacy greedy chaining，也未区分 `BANDED_FAST` / `BANDED_SLOW` 两种延伸模式；因此
  `--default` / `--sensitive` 的延伸比 DIAMOND 更彻底（详见 README「已知差异」）。
- 未实现 `--algo 1`、`--target-indexed`、`--swipe`、`--frameshift`、`--iterate`、`--approx-id`、
  XML/JSON 输出。
- 种子枚举与 DP 均为标量实现，`--sensitive` 及以上比 DIAMOND 慢 4–12 倍。

[Unreleased]: https://github.com/<user>/diamond-rust/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/<user>/diamond-rust/releases/tag/v0.1.0
