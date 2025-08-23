## BioMatics MSA

### 1) Requirements
- Windows 10+ (PowerShell) or Linux/macOS (bash)
- C++17 compiler: g++ or clang++
  - Check: `g++ --version` or `clang++ --version`

### 2) Open the project
Change directory into your project folder (adjust the path to where you cloned/downloaded):
```powershell
cd <your-project-path>  # e.g., C:\path\to\BioMatics
```

### 3) Prepare an example FASTA
We provide `input.fasta.example`. Copy it to `input.fasta` (or place your own FASTA as `input.fasta`):
```powershell
Copy-Item .\input.fasta.example .\input.fasta
```
Note: The program looks for `input.fasta` in the project root. If missing, it will error.

### 4) Allow script execution (if needed)
Temporarily allow PowerShell scripts in the current session:
```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

### 5) Build and run

Windows (PowerShell):
```powershell
./run_msa.ps1
```
The script will:
- Compile: `g++ -g -std=c++17 main_final_syntax_ok.cpp blosum62.cpp -o msa_program.exe`
- Execute: `msa_program.exe`

Linux/macOS (bash):
```bash
chmod +x ./run_msa.sh
./run_msa.sh
```
The script will:
- Compile: `g++ -g -std=c++17 main_final_syntax_ok.cpp blosum62.cpp -o msa_program` (or use clang++)
- Execute: `./msa_program`

### 6) Outputs
- Final alignment: `final_msa.fasta`
- Distance matrix (i,j,EMD): `emd_distance_matrix.csv`
- Detailed log: `msa_debug.log`

### 7) Troubleshooting
- `g++ not found`: Ensure MinGW-w64/MSYS2 is installed and on PATH (Windows). On Linux/macOS, install build tools (e.g., `sudo apt install build-essential` or Xcode Command Line Tools).
- `input.fasta not found`: Copy `input.fasta.example` to `input.fasta` or provide your own.
- Permission errors: Re-run step 4 in the same PowerShell session.


