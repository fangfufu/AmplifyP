# Build the Flet Windows binary and package as a ZIP and Inno Setup installer.
#
# Usage: .\build_windows.ps1

# Source the virtual environment if it exists and is not already active
if ($null -eq $env:VIRTUAL_ENV -and (Test-Path ".venv")) {
    Write-Host "==> Sourcing virtual environment..."
    . .venv\Scripts\Activate.ps1
}

# Disable rich output/animations under CI
if ($env:CI -eq "true") {
    $env:FLET_CLI_NO_RICH_OUTPUT = "1"
}

# Determine the version
$version = $env:VERSION
if ([string]::IsNullOrEmpty($version)) {
    $version = (python -c "import tomllib; print(tomllib.load(open('pyproject.toml', 'rb'))['project']['version'])")
} else {
    # Remove leading 'v' if present (e.g. from git tag)
    $version = $version -replace '^v', ''
}

Write-Host "==> Generating Git SHA..."
python scripts/gen_git_sha.py
if ($LASTEXITCODE -ne 0) { throw "python scripts/gen_git_sha.py failed with exit code $LASTEXITCODE" }

Write-Host "==> Building Flet Windows binary..."
if (Test-Path "build\windows") { Remove-Item -Recurse -Force "build\windows" }
if (Test-Path "build\AmplifyP") { Remove-Item -Recurse -Force "build\AmplifyP" }
flet build windows src -o build/windows --project AmplifyP --yes
if ($LASTEXITCODE -ne 0) { throw "flet build failed with exit code $LASTEXITCODE" }

Write-Host "==> Moving build artefacts..."
Move-Item -Path "build\windows" -Destination "build\AmplifyP"

Write-Host "==> Packaging ZIP archive..."
$zipFile = "amplifyp-windows-$version.zip"
if (Test-Path $zipFile) { Remove-Item -Force $zipFile }
python -c "import shutil; shutil.make_archive('amplifyp-windows-$version', 'zip', 'build', 'AmplifyP')"
if ($LASTEXITCODE -ne 0) { throw "ZIP packaging failed with exit code $LASTEXITCODE" }

# Build the Inno Setup installer if iscc compiler is available
if (Get-Command "iscc" -ErrorAction SilentlyContinue) {
    Write-Host "==> Building Windows installer..."
    iscc amplifyp.iss /DVersion=$version /O.
    if ($LASTEXITCODE -ne 0) { throw "iscc compiler failed with exit code $LASTEXITCODE" }
} else {
    Write-Host "==> Inno Setup (iscc) not found in PATH. Skipping installer build."
    Write-Host "    To install Inno Setup, run: winget install JRSoftware.InnoSetup"
}

Write-Host "==> Build complete: amplifyp-windows-$version.zip"
