# PowerShell script to build and run the MSA program
# Usage:
#   .\run_msa.ps1

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

try {
    # Move to repository root (script directory)
    Set-Location -Path $PSScriptRoot

    # Check for g++
    if (-not (Get-Command g++ -ErrorAction SilentlyContinue)) {
        throw "g++ not found. Please install MinGW-w64 or MSYS2 (with gcc)."
    }

    Write-Host "Compiling (g++ -g -std=c++17)..." -ForegroundColor Cyan
    & g++ -g -std=c++17 -fopenmp main_final_syntax_ok.cpp blosum62.cpp -o msa_program.exe
    if ($LASTEXITCODE -ne 0) {
        throw "Compilation failed with exit code $LASTEXITCODE"
    }
    Write-Host "Compilation successful -> msa_program.exe" -ForegroundColor Green

    Write-Host "Running msa_program.exe..." -ForegroundColor Cyan
    & .\msa_program.exe
    $exitCode = $LASTEXITCODE
    Write-Host "Program exited with code $exitCode"
    exit $exitCode
}
catch {
    Write-Error $_.Exception.Message
    exit 1
}

