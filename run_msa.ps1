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

    $ompFlag = "-fopenmp"
    $cflags = @("-O3","-DNDEBUG","-std=c++17",$ompFlag)
    Write-Host ("Compiling (g++ {0})..." -f ($cflags -join ' ')) -ForegroundColor Cyan
    & g++ @cflags "main_final_syntax_ok.cpp" "blosum62.cpp" "-o" "msa_program.exe"
    if ($LASTEXITCODE -ne 0) {
        throw "Compilation failed with exit code $LASTEXITCODE"
    }
    Write-Host "Compilation successful -> msa_program.exe" -ForegroundColor Green

    # Set reasonable OpenMP defaults if not already set
    if (-not $env:OMP_NUM_THREADS) { $env:OMP_NUM_THREADS = [string][Environment]::ProcessorCount }
    if (-not $env:OMP_PROC_BIND) { $env:OMP_PROC_BIND = "close" }
    if (-not $env:OMP_SCHEDULE) { $env:OMP_SCHEDULE = "static" }

    if (-not (Test-Path "input.fasta")) {
        Write-Warning "input.fasta not found in project root. The program may error if the file is missing."
    }

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

