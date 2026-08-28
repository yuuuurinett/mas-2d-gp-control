param(
    [int]$StartSeed = 1,
    [int]$EndSeed = 50
)

$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$matlab = 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe'
$resultFolder = Join-Path $repo 'result\mc50_all_modes_T20_U0p4_D0p1_M600_structured_R10'
$logFolder = Join-Path $resultFolder 'logs'
New-Item -ItemType Directory -Force -Path $logFolder | Out-Null

for ($seed = $StartSeed; $seed -le $EndSeed; $seed++) {
    $logFile = Join-Path $logFolder ('seed_{0:d3}.log' -f $seed)
    $command = "cd('$($repo.Replace("'","''"))'); " +
        "addpath(genpath(pwd)); " +
        "run_all_modes_mc50_T20(true,$seed,false);"
    Write-Host ("Starting final all-mode seed {0} (requested {1}-{2})." -f $seed, $StartSeed, $EndSeed)
    $oldPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    & $matlab -batch $command *>&1 | Tee-Object -FilePath $logFile
    $matlabExitCode = $LASTEXITCODE
    $ErrorActionPreference = $oldPreference
    if ($matlabExitCode -ne 0) {
        throw "MATLAB failed at seed $seed. See $logFile"
    }
}

$summaryLog = Join-Path $logFolder 'summary.log'
$summaryCommand = "cd('$($repo.Replace("'","''"))'); " +
    "addpath(genpath(pwd)); " +
    "run_all_modes_mc50_T20(false,$($StartSeed):$($EndSeed),true,false);"
$oldPreference = $ErrorActionPreference
$ErrorActionPreference = 'Continue'
& $matlab -batch $summaryCommand *>&1 | Tee-Object -FilePath $summaryLog
$matlabExitCode = $LASTEXITCODE
$ErrorActionPreference = $oldPreference
if ($matlabExitCode -ne 0) {
    throw "MATLAB summary failed. See $summaryLog"
}
Write-Host ("Final all-mode seeds {0}-{1} completed." -f $StartSeed, $EndSeed)
