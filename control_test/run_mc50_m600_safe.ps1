param(
    [int]$StartSeed = 17,
    [int]$EndSeed = 50
)

$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$matlab = 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe'
$logFolder = Join-Path $repo 'result\mc50_all_aggregation_m600_T30_per_agent\logs'
New-Item -ItemType Directory -Force -Path $logFolder | Out-Null

for ($seed = $StartSeed; $seed -le $EndSeed; $seed++) {
    $logFile = Join-Path $logFolder ('seed_{0:d3}.log' -f $seed)
    $matlabCommand = "cd('$($repo.Replace("'","''"))'); " +
        "addpath(genpath(pwd)); " +
        "run_mc10_all_aggregation_methods(true,$seed,false,600,true,50);"

    Write-Host ("Starting M600 seed {0}. Existing MAT files will be skipped." -f $seed)
    & $matlab -batch $matlabCommand *>&1 | Tee-Object -FilePath $logFile
    if ($LASTEXITCODE -ne 0) {
        throw "MATLAB failed at seed $seed. See $logFile"
    }

    # Every seed gets a fresh MATLAB process, releasing all accumulated RAM.
    Start-Sleep -Seconds 3
}

Write-Host 'MC50 M600 run completed.'
Write-Warning 'Summary was intentionally skipped: audit/rerun the mixed M400 common batch first.'
