param(
    [switch]$Configure,
    [switch]$Clean,
    [switch]$Fresh
)

$ErrorActionPreference = 'Stop'
$workspace = $PSScriptRoot
$aliasPath = Join-Path $env:TEMP 'comb-method-lab1-ascii'
$mingwBin = 'C:\msys64\mingw64\bin'
$env:PATH = "$mingwBin;$env:PATH"

if (-not (Test-Path -LiteralPath $aliasPath)) {
    New-Item -ItemType Junction -Path $aliasPath -Target $workspace | Out-Null
}

Push-Location $aliasPath
try {
    if ($Clean) {
        & cmake --build --preset default --target clean
        exit $LASTEXITCODE
    }

    $cache = Join-Path $aliasPath 'build\CMakeCache.txt'
    if ($Configure -or $Fresh -or -not (Test-Path -LiteralPath $cache)) {
        $args = @('--preset', 'default')
        if ($Fresh) {
            $args += '--fresh'
        }
        & cmake @args
        if ($LASTEXITCODE -ne 0) {
            exit $LASTEXITCODE
        }
    }

    if (-not $Configure -or $Fresh) {
        & cmake --build --preset default
        exit $LASTEXITCODE
    }
}
finally {
    Pop-Location
}
