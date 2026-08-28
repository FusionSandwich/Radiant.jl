param(
    [string]$Julia16 = $env:JULIA16_BIN,
    [string]$Julia110 = $env:JULIA110_BIN
)

$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path

if ([string]::IsNullOrWhiteSpace($Julia16)) {
    $Julia16 = 'julia-1.6.7.exe'
}
if ([string]::IsNullOrWhiteSpace($Julia110)) {
    $Julia110 = 'julia-1.10.12.exe'
}

function Invoke-RadiantQualification {
    param(
        [Parameter(Mandatory=$true)][string]$Executable,
        [Parameter(Mandatory=$true)][string]$ExpectedVersion
    )

    $Command = Get-Command $Executable -ErrorAction SilentlyContinue
    if ($null -eq $Command -and -not (Test-Path $Executable)) {
        throw "Missing Julia $ExpectedVersion executable: $Executable"
    }

    $Actual = (& $Executable --startup-file=no -e 'print(VERSION)').Trim()
    if ($LASTEXITCODE -ne 0) {
        throw "Unable to query Julia executable: $Executable"
    }
    if ($Actual -ne $ExpectedVersion) {
        throw "Expected Julia $ExpectedVersion, found $Actual at $Executable"
    }

    Write-Host "Running Radiant qualification with Julia $Actual"
    & $Executable --startup-file=no --project=$Root `
        (Join-Path $Root 'scripts\qualification\run_julia_qualification.jl')
    if ($LASTEXITCODE -ne 0) {
        throw "Radiant qualification failed for Julia $ExpectedVersion"
    }
}

Invoke-RadiantQualification -Executable $Julia16 -ExpectedVersion '1.6.7'
Invoke-RadiantQualification -Executable $Julia110 -ExpectedVersion '1.10.12'
Write-Host 'Pinned Julia qualification matrix passed.'
