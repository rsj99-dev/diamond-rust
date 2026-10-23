# 三方对比：diamond-rust / DIAMOND 2.2.8 / NCBI BLAST+ 的 blastp 结果与耗时。
#
# 数据可以两种方式给：
#   1. -DataDir <目录>：目录下若干 .faa + 一个 目标序列.fa，先合并成一个大库；
#   2. 不给 -DataDir（且自动探测不到本地 测试\蛋白比蛋白）时，直接用 bench\e2e 的
#      端到端样例（7 860 条 / 2 263 996 aa 库，1 条 497 aa 查询）。
#
# 工作目录默认 C:\diamond-cmp：makeblastdb 的 LMDB 卷无法建在非 ASCII 路径下
# （mdb_env_open 报"系统找不到指定的路径"），所以库和输出都放这里。
param(
    [string]$DataDir = "",
    [string]$WorkDir = "C:\diamond-cmp",
    [string]$Diamond = "",
    [string]$BlastBin = "C:\Program Files\NCBI\blast-2.17.0+\bin",
    [int]$Threads = 8,
    [int]$Repeat = 5,
    [double]$Evalue = 0.001
)

$ErrorActionPreference = "SilentlyContinue"
$parent = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$e2e = Join-Path $PSScriptRoot "e2e"
$rust = Join-Path (Split-Path $PSScriptRoot -Parent) "target\release\diamond-rust.exe"

# 本地保留着 DIAMOND 源码树时，二进制通常在 测试\diamond-2.2.8.exe。
if (-not $Diamond) { $Diamond = Join-Path $parent "测试\diamond-2.2.8.exe" }
if (-not $DataDir) {
    $local = Join-Path $parent "测试\蛋白比蛋白"
    if (Test-Path $local) { $DataDir = $local }
}
# 没装到默认位置时退回 PATH 上的 makeblastdb。
if (-not (Test-Path (Join-Path $BlastBin "blastp.exe"))) {
    $onPath = Get-Command makeblastdb.exe -ErrorAction SilentlyContinue
    if ($onPath) { $BlastBin = Split-Path $onPath.Source }
}
$blastp = Join-Path $BlastBin "blastp.exe"
$makeblastdb = Join-Path $BlastBin "makeblastdb.exe"
$presets = @("default", "fast", "sensitive", "more-sensitive", "very-sensitive")
$fields = @("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore")

foreach ($exe in @($rust, $Diamond, $blastp, $makeblastdb)) {
    if (-not (Test-Path $exe)) { throw "missing executable: $exe" }
}

New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
$dbFa = Join-Path $WorkDir "db_prot.fa"
$query = Join-Path $WorkDir "query.fa"

# --- staging: one library file plus one query file -------------------------
$quant = ""
if ($DataDir) { $quant = Join-Path $DataDir "目标序列.fa" }
$srcs = @()
if ($DataDir -and (Test-Path $DataDir)) {
    $srcs = @(Get-ChildItem -LiteralPath $DataDir -Filter *.faa | ForEach-Object { $_.FullName })
}

if ($srcs.Count -gt 0) {
    $lines = New-Object System.Collections.Generic.List[string]
    foreach ($f in $srcs) { $lines.AddRange([System.IO.File]::ReadAllLines($f)) }
    [System.IO.File]::WriteAllLines($dbFa, $lines, (New-Object System.Text.ASCIIEncoding))
    [System.IO.File]::WriteAllLines($query, [System.IO.File]::ReadAllLines($quant), (New-Object System.Text.ASCIIEncoding))
    Write-Host ("library: {0} seqs, {1:N0} aa; query: {2} seqs" -f
        ($lines | Where-Object { $_.StartsWith(">") }).Count,
        (($lines | Where-Object { -not $_.StartsWith(">") }) -join "").Length,
        ([System.IO.File]::ReadAllLines($query) | Where-Object { $_.StartsWith(">") }).Count)
}
else {
    Copy-Item (Join-Path $e2e "db_prot.fa") $dbFa -Force
    Copy-Item (Join-Path $e2e "query.fa") $query -Force
    Write-Host ("library: bench\e2e\db_prot.fa, query: bench\e2e\query.fa")
}

function Get-Median([double[]]$values) {
    $sorted = $values | Sort-Object
    $mid = [int][math]::Floor($sorted.Count / 2)
    if ($sorted.Count % 2 -eq 1) { return $sorted[$mid] }
    return ($sorted[$mid - 1] + $sorted[$mid]) / 2
}

function Invoke-Timed([string]$exe, [string[]]$argv) {
    $times = @()
    for ($i = 0; $i -lt $Repeat; $i++) {
        $times += (Measure-Command { & $exe @argv 2>$null }).TotalSeconds
    }
    return Get-Median $times
}

function Read-Hits([string]$path) {
    $map = [ordered]@{}
    if (-not (Test-Path $path)) { return $map }
    foreach ($line in [System.IO.File]::ReadAllLines($path)) {
        if (-not $line.Trim()) { continue }
        $f = $line -split "`t"
        if ($f.Count -lt 12) { continue }
        $map[[string]$f[1]] = [pscustomobject]@{
            Length = [int]$f[3]; Pid = [double]$f[2]; Bits = [double]$f[11]
            Evalue = $f[10]; QB = $f[6]; QE = $f[7]; SB = $f[8]; SE = $f[9]
        }
    }
    return $map
}

# --- build ----------------------------------------------------------------
$build = [ordered]@{}
$build["diamond-rust"] = Invoke-Timed $rust (@("makedb", "--in", $dbFa, "--db", (Join-Path $WorkDir "db_rust"), "-p", "$Threads", "--quiet"))
$build["diamond"] = Invoke-Timed $Diamond (@("makedb", "--in", $dbFa, "-d", (Join-Path $WorkDir "db_diamond"), "-p", "$Threads", "--quiet"))
$build["blastp"] = Invoke-Timed $makeblastdb (@("-in", $dbFa, "-dbtype", "prot", "-out", (Join-Path $WorkDir "db_blast")))

# --- search ---------------------------------------------------------------
$hits = [ordered]@{}
$time = [ordered]@{}
$base = @("blastp", "--db", (Join-Path $WorkDir "db_rust"), "--query", $query, "-p", "$Threads", "--outfmt", "6") + $fields + @("--quiet")

foreach ($preset in $presets) {
    $flag = @()
    if ($preset -ne "default") { $flag = @("--$preset") }

    $out = Join-Path $WorkDir "rust_$preset.tsv"
    $argv = $base[0..2] + @("-o", $out) + $base[3..($base.Count - 1)] + $flag
    $time["rust/$preset"] = Invoke-Timed $rust $argv
    $hits["rust/$preset"] = Read-Hits $out

    $out = Join-Path $WorkDir "diamond_$preset.tsv"
    $argv = @("blastp", "-d", (Join-Path $WorkDir "db_diamond"), "-q", $query, "-o", $out,
        "-p", "$Threads", "--outfmt", "6") + $fields + @("--quiet") + $flag
    $time["diamond/$preset"] = Invoke-Timed $Diamond $argv
    $hits["diamond/$preset"] = Read-Hits $out
}

foreach ($seg in @("yes", "no")) {
    $out = Join-Path $WorkDir "blast_seg$seg.tsv"
    $argv = @("-query", $query, "-db", (Join-Path $WorkDir "db_blast"), "-out", $out,
        "-outfmt", ("6 " + ($fields -join " ")), "-evalue", "$Evalue",
        "-num_threads", "$Threads", "-max_target_seqs", "25", "-max_hsps", "1", "-seg", $seg)
    $time["blast/$seg"] = Invoke-Timed $blastp $argv
    $hits["blast/$seg"] = Read-Hits $out
}

# --- report ---------------------------------------------------------------
Write-Host ""
Write-Host "makedb (median of $Repeat)"
foreach ($k in $build.Keys) { Write-Host ("  {0,-14} {1,7:N3} s" -f $k, $build[$k]) }

Write-Host ""
Write-Host ("{0,-22} {1,6} {2,9} {3,9}" -f "run", "hits", "seconds", "top bits")
foreach ($k in $hits.Keys) {
    $top = 0.0
    foreach ($v in $hits[$k].Values) { if ($v.Bits -gt $top) { $top = $v.Bits } }
    $sec = "-"
    if ($time.Contains($k)) { $sec = "{0:N3}" -f $time[$k] }
    Write-Host ("{0,-22} {1,6} {2,9} {3,9}" -f $k, $hits[$k].Count, $sec, $top)
}

$reference = $hits["blast/no"]
if ($reference -and $reference.Count -gt 0) {
    Write-Host ""
    Write-Host "targets recovered vs blastp -seg no"
    foreach ($k in $hits.Keys) {
        $shared = @($hits[$k].Keys | Where-Object { $reference.Contains($_) }).Count
        Write-Host ("  {0,-22} {1,2} / {2,2}" -f $k, $shared, $reference.Count)
    }
}

# Per-target view: length / identity / bits for the three engines.
Write-Host ""
Write-Host ("{0,-18} {1,-20} {2,-20} {3,-20} {4}" -f "subject", "rust/sensitive", "diamond/sens", "diamond/more", "blastp(seg no)")
$runs = @($hits["rust/sensitive"], $hits["diamond/sensitive"], $hits["diamond/more-sensitive"], $hits["blast/no"])
$subjects = @()
foreach ($m in $runs) { $subjects += @($m.Keys) }
foreach ($s in ($subjects | Sort-Object -Unique)) {
    $cells = @()
    foreach ($m in $runs) {
        if ($m.Contains($s)) { $v = $m[$s]; $cells += ("{0}L {1}% {2}" -f $v.Length, $v.Pid, $v.Bits) }
        else { $cells += "--" }
    }
    Write-Host ("{0,-18} {1,-20} {2,-20} {3,-20} {4}" -f $s, $cells[0], $cells[1], $cells[2], $cells[3])
}
