@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\run_examples.ps1" %*
set "HOSC_RESULT=%ERRORLEVEL%"
echo.
if not "%HOSC_RESULT%"=="0" echo Examples did not pass. Read the error above.
pause
exit /b %HOSC_RESULT%
