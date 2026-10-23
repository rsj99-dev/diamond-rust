# Compare diamond-rust with the reference DIAMOND binary.
#
# Builds a database with both programs, runs blastp with the same sensitivity
# presets, prints wall-clock timings and diffs the tabular output line by line.
#
# Usage:
#   .\compare.ps1 -Diamond C:\tools\diamond.exe
#   .\compare.ps1 -Diamond <diamond.exe> -Scale 4 -Repeat 3
#
# -Diamond is the reference DIAMOND 2.2.8 binary; everything else defaults to
# bench\e2e and target\release.
#
# -Scale N repeats the database N times (IDs get a `_copyN` suffix) to measure
# how both programs scale to a larger reference.

param(
    [string]$Diamond,
    [string]$RustBin = "",
    [string]$Fasta = "",
    [string]$Query = "",
    [string[]]$Presets = @("default", "fast", "sensitive", "more-sensitive", "very-sensitive"),
    [int]$Threads = 8,
    [int]$Scale = 1,
    [int]$Repeat = 2
)

$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
if ($Fasta -eq "") { $Fasta = Join-Path $here "e2e\db_prot.fa" }
if ($Query -eq "") { $Query = Join-Path $here "e2e\query.fa" }
if ($RustBin -eq "") { $RustBin = Join-Path (Split-Path -Parent $here) "target\release\diamond-rust.exe" }
if (-not (Test-Path $RustBin)) { throw "diamond-rust not found at $RustBin, run .\build.ps1 first" }

$work = Join-Path $here "e2e\run"
New-Item -ItemType Directory -Force -Path $work | Out-Null

$fmt = @("6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore")

function Expand-Fasta([string]$inPath, [string]$outPath, [int]$copies) {
    $lines = [System.IO.File]::ReadAllLines($inPath)
    $writer = New-Object System.IO.StreamWriter($outPath)
    try {
        for ($c = 1; $c -le $copies; $c++) {
            foreach ($line in $lines) {
                if ($line.StartsWith(">")) { $writer.WriteLine($line + "_copy$c") }
                else { $writer.WriteLine($line) }
            }
        }
    }
    finally { $writer.Close() }
}

function Invoke-Timed([string]$exe, [string[]]$arguments, [int]$repeats) {
    # Native stderr is reported as a PowerShell error; quiet mode plus a local
    # preference keeps it from aborting the script.
    $ErrorActionPreference = "SilentlyContinue"
    $best = [double]::MaxValue
    for ($i = 0; $i -lt $repeats; $i++) {
        $t = (Measure-Command { & $exe @arguments 2>$null | Out-Null }).TotalSeconds
        if ($t -lt $best) { $best = $t }
    }
    return $best
}

# Tabular hits keyed by `qseqid|sseqid`; --max-hsps 1 means one row per key.
function Read-Hits([string]$path) {
    $map = @{}
    foreach ($line in [System.IO.File]::ReadAllLines($path)) {
        if ($line.Trim() -eq "") { continue }
        $f = $line -split "`t"
        if ($f.Count -lt 12) { continue }
        $map["$($f[0])|$($f[1])"] = [pscustomobject]@{
            Line   = $line
            QBegin = $f[6]; QEnd = $f[7]; SBegin = $f[8]; SEnd = $f[9]
        }
    }
    return $map
}

# --- database -----------------------------------------------------------------
$dbFasta = $Fasta
if ($Scale -gt 1) {
    $dbFasta = Join-Path $work "db_x$Scale.fa"
    Write-Host "Expanding $Fasta x$Scale -> $dbFasta"
    Expand-Fasta $Fasta $dbFasta $Scale
}

$refDb = Join-Path $work "ref"
$rustDb = Join-Path $work "rust"

$rustMakedb = @("makedb", "--in", $dbFasta, "--db", $rustDb, "-p", $Threads, "--quiet")
$tRustDb = Invoke-Timed $RustBin $rustMakedb 1
Write-Host ("makedb  diamond-rust  {0,7:N3}s" -f $tRustDb)

$tRefDb = $null
if ($Diamond) {
    if (-not (Test-Path $Diamond)) { throw "diamond not found at $Diamond" }
    $refMakedb = @("makedb", "--in", $dbFasta, "-d", $refDb, "-p", $Threads, "--quiet")
    $tRefDb = Invoke-Timed $Diamond $refMakedb 1
    Write-Host ("makedb  diamond       {0,7:N3}s" -f $tRefDb)
}

# --- searches -----------------------------------------------------------------
$rows = @()
foreach ($preset in $Presets) {
    $flag = $null
    if ($preset -ne "default") { $flag = "--$preset" }

    $rustOut = Join-Path $work "rust_$preset.tsv"
    $rustArgs = @("blastp", "--db", $rustDb, "--query", $Query, "--out", $rustOut,
        "-p", $Threads, "--outfmt") + $fmt + @("--quiet")
    if ($flag) { $rustArgs += $flag }
    $tRust = Invoke-Timed $RustBin $rustArgs $Repeat

    $diamondRows = 0
    $rustRows = 0
    $shared = 0
    $coords = 0
    $exact = 0
    $tRef = $null
    if ($Diamond) {
        $refOut = Join-Path $work "ref_$preset.tsv"
        $refArgs = @("blastp", "-d", $refDb, "-q", $Query, "-o", $refOut,
            "-p", $Threads, "--outfmt") + $fmt + @("--quiet")
        if ($flag) { $refArgs += $flag }
        $tRef = Invoke-Timed $Diamond $refArgs $Repeat
        $refHits = Read-Hits $refOut
        $rustHits = Read-Hits $rustOut
        $diamondRows = $refHits.Count
        $rustRows = $rustHits.Count
        $common = @($refHits.Keys | Where-Object { $rustHits.ContainsKey($_) })
        $shared = $common.Count
        $coords = @($common | Where-Object {
                $d = $refHits[$_]; $r = $rustHits[$_]
                $d.QBegin -eq $r.QBegin -and $d.QEnd -eq $r.QEnd -and
                $d.SBegin -eq $r.SBegin -and $d.SEnd -eq $r.SEnd
            }).Count
        $exact = @($common | Where-Object { $refHits[$_].Line -eq $rustHits[$_].Line }).Count
        $onlyRef = @($refHits.Keys | Where-Object { -not $rustHits.ContainsKey($_) })
        $onlyRust = @($rustHits.Keys | Where-Object { -not $refHits.ContainsKey($_) })
        foreach ($k in @($onlyRef | Select-Object -First 3)) {
            Write-Host ("  diamond only [{0}]: {1}" -f $preset, ($refHits[$k].Line -replace "`t", " | "))
        }
        foreach ($k in @($onlyRust | Select-Object -First 3)) {
            Write-Host ("  rust only [{0}]: {1}" -f $preset, ($rustHits[$k].Line -replace "`t", " | "))
        }
    }
    else {
        $rustRows = (Read-Hits $rustOut).Count
        $diamondRows = $rustRows
    }

    $dSec = $null
    if ($null -ne $tRef) { $dSec = [math]::Round($tRef, 3) }
    $rows += [pscustomobject]@{
        Preset      = $preset
        DiamondHits = $diamondRows
        RustHits    = $rustRows
        Shared      = $shared
        Coords      = $coords
        Exact       = $exact
        DiamondSec  = $dSec
        RustSec     = [math]::Round($tRust, 3)
    }
}

Write-Host ""
Write-Host "preset            hitsD hitsR  shared  coords   exact   t_diamond   t_rust"
foreach ($r in $rows) {
    $td = "-"
    if ($null -ne $r.DiamondSec) { $td = "{0,8:N3}s" -f $r.DiamondSec }
    Write-Host ("{0,-16} {1,5} {2,5} {3,7} {4,7} {5,7}   {6}  {7,7:N3}s" -f `
            $r.Preset, $r.DiamondHits, $r.RustHits, $r.Shared, $r.Coords, $r.Exact, $td, $r.RustSec)
}
Write-Host ""
Write-Host "Outputs and databases: $work"
if ($Scale -gt 1) { Write-Host "Scaled FASTA: $dbFasta" }
