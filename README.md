# diamond-rust

[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)
[![Rust 1.75+](https://img.shields.io/badge/rust-1.75%2B-orange.svg)](Cargo.toml)
![Platforms](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)
<!-- 建仓后把 <user> 换成 GitHub 用户名，再取消下一行的注释。
[![CI](https://github.com/<user>/diamond-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/<user>/diamond-rust/actions/workflows/ci.yml)
-->

**diamond-rust** 是用 Rust 重写的 [DIAMOND](https://github.com/bbuchfink/diamond)（Double Index
Aligner）。命令行、灵敏度预设与 `--outfmt 6` 输出对齐 DIAMOND 2.2.8，算法流水线按 DIAMOND 的
C++ 源码（`diamond-2.2.8/src`，仅本地参考、不随仓库分发）逐段移植：

> makedb → 双索引种子检索（约简字母表 × spaced seed）→ ungapped X-drop → 对角线分段
> → 带状仿射 Smith-Waterman + traceback → ALP E-value/bitscore → 排序/culling → outfmt 6

在 2.26 Maa 的蛋白库、497 aa 查询、8 线程下：`--fast` 与 `--default` 的耗时接近 DIAMOND
（0.23 s vs 0.20 s，0.45 s vs 0.25 s），`--sensitive` 及以上慢 4–12 倍（瓶颈是标量种子哈希，
DIAMOND 用 SIMD）；`--sensitive`/`--more-sensitive`/`--very-sensitive` 的命中目标集合与
DIAMOND 完全一致，且与 NCBI BLAST+ 的 7 个远缘同源全部对得上。详见[保真度与性能](#保真度与性能)。

## 支持的命令

| 命令 | 查询 | 数据库 |
| --- | --- | --- |
| `makedb` | FASTA | 写出 `*.dmnd` |
| `blastp` | 蛋白 | 蛋白 |
| `blastx` | 核酸（六帧翻译） | 蛋白 |
| `dbinfo` | — | 打印库统计 |

数据库是本工具自己的二进制格式（magic `DRDMND01`），**不兼容 DIAMOND 的 `.dmnd`**：
换工具就要重新 `makedb`。格式只存序列与 ID，种子索引在搜索时于内存中构建，与 DIAMOND
2.x 的 `--algo 0` 一致。

## 安装与构建

需要 Rust 1.75 及以上，只依赖 `rayon`。

```bash
cargo install --path .         # 安装到 ~/.cargo/bin
exe="target/release/diamond-rust"; cargo build --release   # 或只构建：$exe / diamond-rust.exe
```

Windows 上默认的 `x86_64-pc-windows-msvc` 需要 Visual Studio 生成工具里的 C++ 桌面开发组件
（提供 `link.exe`）。没有该组件时，本仓库的 [`build.ps1`](build.ps1) 会自动回退到
`stable-x86_64-pc-windows-gnu` + Rust 自带的 `rust-lld`，无需另装 MinGW/gcc：

```powershell
.\build.ps1              # release 构建
.\build.ps1 -Test        # 构建 + 单元测试（57 个）
.\build.ps1 -Debug       # debug 构建
```

等价的手工命令：

```powershell
rustup toolchain install stable-x86_64-pc-windows-gnu
$env:CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER = "rust-lld"
$env:RUSTFLAGS = "-Clink-self-contained=yes"
cargo +stable-x86_64-pc-windows-gnu build --release
```

`-Clink-self-contained=yes` 是稳定选项，链接产物不依赖 `libgcc_s_seh-1.dll` 等外部 DLL。

## 快速开始

```bash
# 1. 建库（一次性）
diamond-rust makedb --in db_prot.fa --db db_prot

# 2. 搜索（不写 --out 时输出到标准输出）
diamond-rust blastp --db db_prot.dmnd --query query.fa --outfmt 6 --out hits.tsv

# 3. 看库信息
diamond-rust dbinfo --db db_prot.dmnd
```

`blastx` 的用法与 `blastp` 相同，查询换成核酸 FASTA；输出坐标会换算回核酸序列（`qframe`
为 1–3 / -1–-3），`qcovhsp` 按核酸长度计算，与 DIAMOND 一致。

## 常用参数

| 参数 | 说明 |
| --- | --- |
| `-e, --evalue <VALUE>` | 报告阈值，默认 `0.001` |
| `-k, --max-target-seqs <N>` | 每个查询最多保留的目标数，默认 25 |
| `--max-hsps <N>` | 每个目标最多保留的 HSP 数，默认 1 |
| `-p, --threads <N>` | 线程数，默认为全部核心 |
| `--outfmt 6 [field ...]` | 见表下字段列表 |
| `--matrix <NAME>` | 打分矩阵，默认 `BLOSUM62`（内置 BLOSUM/PAM 系列） |
| `--gapopen` / `--gapextend` | 空位罚分，默认随矩阵（BLOSUM62 为 11/1） |
| `--seed-cut <VALUE>` | 种子复杂度阈值，默认随灵敏度预设（0.8–1.0） |
| `--index-chunks <N>` | 一次索引的库分块数，默认随预设（`very-sensitive`/`ultra-sensitive` 为 1） |
| `--fast` … `--ultra-sensitive` | 灵敏度预设，与 DIAMOND 同名同义 |

灵敏度预设的 shape 集、约简字母表、seed 复杂度阈值、band 宽度、`index_chunks` 均取自
`search/setup.cpp` 与 `align/gapped_score.cpp`，不是经验值。例如 `--very-sensitive` 用 14×7
的 16 个 shape、`diamond9` 约简、整库单块索引；其余预设用 `murphy10`，只对同聚物种子做屏蔽。

`--outfmt` 支持的字段：`qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue
bitscore score length pident nident mismatch positive gapopen gaps ppos qframe qcovhsp scovhsp
cigar btop qtitle stitle full_sseq qseq_translated qstrand`。数字格式对齐 DIAMOND：
E-value 用 `%.2e`，bitscore/pident/ppos/qcovhsp/scovhsp 在 ≥100 时取整、<100 时保留一位小数。

DIAMOND 的其余开关（`--algo`、`--block-size`、`--comp-based-stats`、`--masking` 等）会被接受并
忽略，因此同一串 DIAMOND 命令行可以直接跑；但依赖这些开关的功能并不生效，见[已知差异](#已知差异)。

## 代码结构

| 文件 | 职责 |
| --- | --- |
| `src/cli.rs` | 参数解析（`makedb`/`blastp`/`blastx` 的 DIAMOND 兼容长选项） |
| `src/io.rs` | FASTA 读取与核酸/蛋白识别 |
| `src/alphabet.rs` | 20 字母编码（NCBI 顺序）、`X` 与非法字符处理 |
| `src/reduction.rs` | 约简字母表：`murphy10`、`diamond9`/`diamond14`/`diamond18` |
| `src/shape.rs` | spaced seed shape 的解析与打包（按约简字母表定宽） |
| `src/seed.rs` | 种子枚举、低复杂度过滤、query/ref 双索引归并 |
| `src/ungapped.rs` | 对角线 X-drop 延伸，输出 `Segment` |
| `src/gapped.rs` | 带状仿射 Smith-Waterman + traceback（`extend`/`extend_span`） |
| `src/matrix.rs`、`matrix_data.rs` | 打分矩阵、Karlin-Altschul 参数、ALP E-value、bit score |
| `src/hsp.rs` | HSP 结构、坐标与统计、CIGAR/BTOP |
| `src/output.rs` | `--outfmt 6` 字段渲染与数字格式 |
| `src/search.rs` | 流水线编排：shape/约简/seed-cut 预设、band 合并、culling |
| `src/translate.rs` | `blastx` 六帧翻译、`--min-orf` 与简并密码子处理 |
| `src/db.rs` | `.dmnd` 二进制库格式（magic `DRDMND01`，v1） |
| `src/workflow.rs` | 子命令分发与端到端流程 |

## 保真度与性能

测试环境：Windows 11 / 8 线程；库 `bench/e2e/db_prot.fa`（7 860 条、2 263 996 aa），查询
`bench/e2e/query.fa`（1 条、497 aa）；对比对象是 DIAMOND 2.2.8 的官方二进制与 NCBI BLAST+
2.17.0 的 `blastp`（`-evalue 1e-3 -max_target_seqs 25 -max_hsps 1`）。计时为 5 次取中位数，
DIAMOND 一侧用 `--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart
send evalue bitscore --quiet`，BLAST 用同字段的 `-outfmt "6 ..."`。

| 预设 | DIAMOND 命中 | diamond-rust 命中 | blastp 命中 | DIAMOND 耗时 | diamond-rust 耗时 |
| --- | --- | --- | --- | --- | --- |
| `--default` | 3 | 4 | 7 | 0.25 s | 0.45 s |
| `--fast` | 1 | 1 | 7 | 0.20 s | 0.23 s |
| `--sensitive` | 7 | 7 | 7 | 0.98 s | 4.1 s |
| `--more-sensitive` | 7 | 7 | 7 | 0.95 s | 4.2 s |
| `--very-sensitive` | 7 | 7 | 7 | 0.49 s | 5.8 s |

`blastp` 没有灵敏度档位，单条 497 aa 查询耗时 2.0 s（`-seg no`）/ 2.8 s（`-seg yes`）；
建库 2.9 s、4.5 MB，而 DIAMOND 与 diamond-rust 的 `makedb` 都在 0.06 s 内（2.9 MB）。
表里的耗时为端到端墙钟时间（含进程启动与输出写入），由 `bench/compare3.ps1` 测得。

`--sensitive` 及以上三方的 7 个目标完全相同，本实现各条 HSP 的列数与 DIAMOND
`--more-sensitive` 相差不超过 2 列（`WP_048121137.1` 两边都是 612 列）。同一段比对三方的
bitscore 一致（例：482 列的对齐都是 321），差异来自 HSP 边界而非打分公式。

复现：

```powershell
.\bench\compare.ps1 -Diamond C:\tools\diamond.exe -Threads 8 -Repeat 3
.\bench\compare.ps1 -Diamond C:\tools\diamond.exe -Scale 4    # 库放 4 倍，看扩展性
.\bench\compare3.ps1 -Diamond C:\tools\diamond.exe -Repeat 5 # 加入 NCBI BLAST+ 的三方对比
```

两个脚本默认都读 `bench/e2e`；`compare3.ps1` 另外接受 `-DataDir <目录>`（一组 `.faa`
加一份 `目标序列.fa`，会先合并成一个大库）。它默认在 `C:\diamond-cmp` 下建库并输出：
`makeblastdb` 的 LMDB 卷必须在纯 ASCII 路径下创建，否则报 `mdb_env_open`。BLAST+ 默认到
`C:\Program Files\NCBI\blast-2.17.0+\bin` 找，也可以用 `-BlastBin` 指定。

## 已知差异

1. **未实现 DIAMOND 的 legacy greedy chaining，也没有区分 `BANDED_FAST` /
   `BANDED_SLOW` 两种延伸模式**（`legacy/chaining/greedy_align.cpp`、
   `align/extend.cpp:102-115`）。替代做法是把同一 target 上所有 `[diag-band, diag+band]`
   重叠的 ungapped 段合并成一次宽带 DP（`align/gapped_score.cpp`），锚点取组内最高分段
   的起点。DIAMOND 在 `--default`～`--sensitive` 用窄带快速模式（497 aa 查询的带宽 64），
   `--more-sensitive` 及以上才切到宽带（120）；本实现每一档都是合并后的宽带完整仿射 SW，
   等效于后者。后果是 `--default` 多报 `WP_048121132.1`（26.2%、83.2 bits），`--sensitive`
   下 `WP_048121137.1` 直接给出 612 列、91.7 bits，而 DIAMOND `--sensitive` 只给 144 列、
   78.6 bits。已验证 `diamond blastp --default --band 120` 能复现 612 列，说明这一处差异
   只来自带宽；`--masking 0`、`--freq-sd 999` 都不改变 DIAMOND 的结果。
2. **HSP 划分不同**：对 `WP_048121137.1`，BLAST 给出两个 HSP（236 列 raw 205 + 192 列
   raw 118，query 4-211 与 304-478），DIAMOND 的宽带模式与本实现把它们缝成一条 612 列
   比对（query 4-478）。长比对不代表找得更好，而是拼接策略不同。
3. **分数不完全相同**：锚点选择与 DP 起点的细微差别让同一 HSP 的 raw score 相差 0–5，
   bitscore 相差 0–2.4，E-value 因此在量级内小偏。
4. **未实现的开关**：`--algo 1`（query-indexed）、`--target-indexed`、`--swipe`、
   `--frameshift`、`--iterate`、`--approx-id`/`--id`、`--query-cover`/`--subject-cover`、
   `--toppercent`/`--top`、XML/JSON 输出、`--masking` 与 `--soft-masking`。
5. **性能上限**：种子枚举是标量循环，DIAMOND 用 SIMD（AVX2/AVX-512）；库越大差距越明显。
   本实现的 DP 也是标量带状 SW，没有 SIMD swipe。

## 测试

```powershell
.\build.ps1 -Test        # 57 个单测：字母表、约简、shape、seed、X-drop、SW、矩阵、输出、CLI
cargo test --release     # 或直接跑
```

单测覆盖的关键不变量：seed shape 编码与 DIAMOND 预设一致、`murphy10`/`diamond9` 分组、
seed 复杂度阈值与 `seed_cut * ln2 * weight` 的换算、band 宽度分档、ALP E-value 与
Karlin-Altschul 公式、`format_double`/`%.2e` 的输出格式、`blastx` 反向框的坐标映射。

## 仓库结构

```
diamond-rust/
├── src/                    # 实现，逐文件说明见上表
├── bench/
│   ├── e2e/                # db_prot.fa（7 860 条）、query.fa、subject.fa
│   ├── compare.ps1         # 与 DIAMOND 对拍：shared / coords / exact 三层统计
│   └── compare3.ps1        # 加入 NCBI BLAST+ 的三方对比
├── scripts/gen_matrices.py # 从 DIAMOND 源码生成 src/matrix_data.rs
├── .github/                # CI 工作流与 issue / PR 模板
├── build.ps1               # Windows 构建，无 link.exe 时回退到 GNU + rust-lld
└── Cargo.toml
```

不随仓库分发的内容：`diamond-2.2.8/`（C++ 参考源码，仅本地保留）、`target/`、
`bench/e2e/run/` 与建库产物 `*.dmnd`。前两者在 `.gitignore` 里，`cargo package` 另外用
`exclude` 把它们和 `bench/e2e` 的大文件排除在 `.crate` 之外。

## 贡献

流程与代码约定见 [CONTRIBUTING.md](CONTRIBUTING.md)。提 PR 前请确认这三条命令全绿，
并在描述里附上与 DIAMOND 的对拍结果：

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --release
```

版本变更见 [CHANGELOG.md](CHANGELOG.md)，安全问题按 [SECURITY.md](SECURITY.md) 私密报告。

## 许可证

[GPL-3.0-or-later](LICENSE)。DIAMOND 本身以 GPL-3.0 分发，本项目的灵敏度预设、打分矩阵
与 ALP 参数又直接移植自其源码，因此沿用同一许可证。

## 致谢

- [DIAMOND](https://github.com/bbuchfink/diamond)（Benjamin Buchfink 等）——算法流程、
  灵敏度预设与输出格式的来源，本实现的可对拍性完全建立在它的源码之上；
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  ——作为第三方基准，用来确认远缘同源没有漏报；
- [rayon](https://github.com/rayon-rs/rayon)——查询级数据并行。

引用信息见 [CITATION.cff](CITATION.cff)。
