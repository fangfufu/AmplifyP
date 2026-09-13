# Build the Flet Windows binary and package as a ZIP and Inno Setup installer.
#
# Usage: .\build_windows.ps1

# Source the virtual environment if it exists and is not already active
if ($null -eq $env:VIRTUAL_ENV -and (Test-Path ".venv")) {
    Write-Host "==> Sourcing virtual environment..."
    . .venv\Scripts\Activate.ps1
}

# Disable rich output and skip Flutter Doctor Android toolchain checks
$env:FLET_CLI_NO_RICH_OUTPUT = "1"
$env:FLET_CLI_SKIP_FLUTTER_DOCTOR = "1"

# Disable Android toolchain requirement in Flutter config if flutter command is available
if (Get-Command "flutter" -ErrorAction SilentlyContinue) {
    flutter config --no-enable-android | Out-Null
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
flet build windows src -o build/windows --project AmplifyP --build-version $version --yes
if ($LASTEXITCODE -ne 0) { throw "flet build failed with exit code $LASTEXITCODE" }

Write-Host "==> Moving build artefacts..."
Move-Item -Path "build\windows" -Destination "build\AmplifyP"

# Clean up temporary files
if (Test-Path "src\amplifyp\gui\git_sha.py") { Remove-Item -Force "src\amplifyp\gui\git_sha.py" }
if (Test-Path "src\build") { Remove-Item -Recurse -Force "src\build" }
if (Test-Path "src\dist") { Remove-Item -Recurse -Force "src\dist" }

Write-Host "==> Packaging ZIP archive..."
$zipFile = "amplifyp-windows-$version.zip"
if (Test-Path $zipFile) { Remove-Item -Force $zipFile }
python -c "import shutil; shutil.make_archive('amplifyp-windows-$version', 'zip', 'build', 'AmplifyP')"
if ($LASTEXITCODE -ne 0) { throw "ZIP packaging failed with exit code $LASTEXITCODE" }

# Build the Inno Setup installer if iscc compiler is available
$isccPath = $null
if (Get-Command "iscc" -ErrorAction SilentlyContinue) {
    $isccPath = "iscc"
} else {
    $knownIsccPaths = @(
        "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
        "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
        "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
    )
    foreach ($path in $knownIsccPaths) {
        if (Test-Path $path) {
            $isccPath = $path
            break
        }
    }
}

if ($isccPath) {
    Write-Host "==> Building Windows installer with $isccPath..."
    & $isccPath amplifyp.iss /DVersion=$version /O.
    if ($LASTEXITCODE -ne 0) { throw "iscc compiler failed with exit code $LASTEXITCODE" }
    Write-Host "==> Build complete: $zipFile and amplifyp-windows-setup-$version.exe"
} else {
    Write-Host "==> Inno Setup (iscc) not found in PATH or standard directories. Skipping installer build."
    Write-Host "    To install Inno Setup, run: winget install JRSoftware.InnoSetup"
    Write-Host "==> Build complete: $zipFile"
}
