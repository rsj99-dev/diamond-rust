# 贡献指南

`diamond-rust` 的目标不是"再写一个比对器"，而是**在不看结果的情况下也能对上台 DIAMOND 的输出**。
所以下面这些约定比一般项目更严格一点。

## 环境

- Rust **1.75+**（`Cargo.toml` 里的 `rust-version`）。
- Windows 上没有 MSVC 的 `link.exe` 时，用 GNU 工具链 + `rust-lld`：

  ```powershell
  rustup toolchain install stable-x86_64-pc-windows-gnu
  $env:CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER = "rust-lld"
  $env:RUSTFLAGS = "-Clink-self-contained=yes"
  cargo +stable-x86_64-pc-windows-gnu build --release
  ```

  仓库自带的 `build.ps1` 会自动做这个回退：

  ```powershell
  .\build.ps1          # release 构建
  .\build.ps1 -Test    # 构建并跑单测
  ```

## 提交前必须通过

```powershell
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --release
```

CI（`.github/workflows/ci.yml`）会在 Linux / macOS / Windows 上重跑这三条，任何一条红都不会合入。

## 改算法的约定

1. **每一处启发式都要能指回 C++ 源码**。注释里写清文件与函数，例如：

   ```rust
   // `setup.cpp`: `cfg.seed_complexity_cut = seed_cut * ln(2) * shapes[0].weight`.
   ```

   没有出处的"经验参数"一律视为 bug。

2. **先对拍，再改**。基准数据在 `bench/e2e/`（7 860 条蛋白、497 aa 查询），两个脚本：

   ```powershell
   .\bench\compare.ps1 -Diamond ..\测试\diamond-2.2.8.exe -Repeat 3   # 与 DIAMOND 对比
   .\bench\compare3.ps1 -Repeat 5                                     # 再加 NCBI BLAST+
   ```

   PR 描述里附上改前 / 改后的命中数、坐标一致数和耗时。

3. **不要为了"多报一条"而放宽阈值**。DIAMOND 的灵敏度分档是有意为之，改动要么对齐预设，
   要么新增一个显式开关。

4. **数值常量原样保留**。`src/matrix_data.rs` 和 `src/alp.rs` 里的系数由
   `scripts/gen_matrices.py` 从 DIAMOND 的 `src/stats/matrices/*.h` 生成，不要手改。

## 代码风格

- `rustfmt` 默认配置，不另加 `rustfmt.toml`。
- 注释解释**为什么**，不复述代码在做什么。
- 错误信息面向使用者：`error: database not found: <path>`，而不是 `unwrap()`。
- 新增 CLI 开关时，同时更新 `src/workflow.rs` 的 usage 文本和 `README` 的参数表。

## 提交信息

用祈使句、首字母大写、不加句号：

```
Fix band merging across touching diagonals
Add --index-chunks override
```

一次提交只做一件事。格式化改动单独提交。

## 报告问题

- 功能异常 / 结果不一致：用 Bug report 模板，附 `diamond-rust version`、命令行、输入大小和
  期望输出（最好有 DIAMOND 的对照结果）。
- 安全问题：见 [SECURITY.md](SECURITY.md)，不要开公开 issue。
