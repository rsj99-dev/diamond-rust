# Build diamond-rust.
#
# On Windows the default `x86_64-pc-windows-msvc` target needs `link.exe` from
# the Visual Studio C++ build tools. When that is missing, this script falls
# back to the GNU toolchain and links with Rust's bundled `rust-lld`, so no
# MinGW installation is required.
#
# Usage:
#   .\build.ps1              # release build
#   .\build.ps1 -Test        # release build + unit tests
#   .\build.ps1 -Debug       # debug build

param(
    [switch]$Test,
    [switch]$Debug
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $root
try {
    # `$PROFILE` is a PowerShell automatic variable, so the build flags need
    # a different name; an empty array splats away cleanly.
    $buildArgs = @()
    if (-not $Debug) { $buildArgs += "--release" }
    # `cargo +<toolchain>` must precede the subcommand; keep it as one word so
    # PowerShell does not hand the `+` to the shell as a separate argument.
    $toolchain = $null

    if (-not (Get-Command link.exe -ErrorAction SilentlyContinue)) {
        $gnu = "stable-x86_64-pc-windows-gnu"
        $installed = & rustup toolchain list 2>$null
        if (-not ($installed -match [regex]::Escape($gnu))) {
            Write-Host "Installing $gnu ..."
            & rustup toolchain install $gnu
        }
        $env:CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER = "rust-lld"
        $env:RUSTFLAGS = "-Clink-self-contained=yes"
        $toolchain = "+$gnu"
        Write-Host "link.exe not found: building with $gnu and rust-lld"
    }

    if ($toolchain) {
        & cargo $toolchain build @buildArgs
    } else {
        & cargo build @buildArgs
    }
    if ($LASTEXITCODE -ne 0) { throw "cargo build failed" }

    if ($Test) {
        if ($toolchain) {
            & cargo $toolchain test @buildArgs
        } else {
            & cargo test @buildArgs
        }
        if ($LASTEXITCODE -ne 0) { throw "cargo test failed" }
    }

    $exe = if ($Debug) { "target\debug\diamond-rust.exe" } else { "target\release\diamond-rust.exe" }
    Write-Host "Built $root\$exe"
}
finally {
    Pop-Location
}
