# Export existing author MATLAB .mat workspaces without recomputing GE.
# Run from any Windows PowerShell directory:
#   powershell -ExecutionPolicy Bypass -File .\scripts\export_author_reference.ps1
param(
    [string]$MatDirectory = 'C:\github\HANKSOME',
    [string]$MatlabExecutable = 'matlab'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$exporterDir = Join-Path $projectRoot 'matlab'
$exporter = Join-Path $exporterDir 'export_author_reference_fixtures.m'
$outputDir = Join-Path (Join-Path $projectRoot 'results') 'author_reference'

if (-not (Test-Path -LiteralPath $exporter -PathType Leaf)) {
    throw "MATLAB exporter missing: $exporter"
}
$inputNames = @('exp_flex.mat', 'exp_flex_lev16.mat', 'exp_both_41_lev16.mat')
foreach ($name in $inputNames) {
    $source = Join-Path $MatDirectory $name
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
        throw "Original author solution missing: $source"
    }
}

$null = Get-Command $MatlabExecutable -ErrorAction Stop
New-Item -ItemType Directory -Path $outputDir -Force | Out-Null

# Use resolved ABSOLUTE paths, never MATLAB pwd/startup working directory.
# Forward slashes are accepted by MATLAB on Windows. Double embedded single
# quotes for MATLAB string literals.
function To-MatlabLiteral([string]$path) {
    $full = [System.IO.Path]::GetFullPath($path).Replace('\', '/')
    return "'" + $full.Replace("'", "''") + "'"
}
$matlabCode = 'addpath(' + (To-MatlabLiteral $exporterDir) + '); ' +
              'export_author_reference_fixtures(' +
              (To-MatlabLiteral $MatDirectory) + ',' +
              (To-MatlabLiteral $outputDir) + ');'
Write-Host 'Exporting ORIGINAL author MATLAB saved solutions (no fsolve)...'
& $MatlabExecutable -batch $matlabCode
if ($LASTEXITCODE -ne 0) {
    throw "MATLAB fixture export failed (exit code $LASTEXITCODE)."
}

$expected = @('original_boom_candidate.csv',
              'shock_pre_distribution.csv',
              'shock_inherited_metadata.csv')
foreach ($label in @('unleveraged', 'leveraged', 'fixed')) {
    foreach ($suffix in @('candidate', 'source_residual', 'source_terms')) {
        $expected += "$label`_$suffix.csv"
    }
}
foreach ($filename in $expected) {
    $f = Join-Path $outputDir $filename
    if (-not (Test-Path -LiteralPath $f -PathType Leaf)) {
        throw "MATLAB returned success but required CSV is missing: $f"
    }
    if ((Get-Item -LiteralPath $f).Length -eq 0) {
        throw "MATLAB fixture is empty: $f"
    }
}
Write-Host "Export OK: $($expected.Count) required nonempty CSVs in $outputDir"
Write-Host 'Next: julia --project=. scripts/evaluate_author_matlab_solution.jl results/author_reference'
