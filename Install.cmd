@echo off
setlocal
echo HOSC plugin installer for 32-bit SGeMS 2.0
echo Save your project and close the target SGeMS before continuing.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\install_plugin.ps1" %*
set "HOSC_INSTALL_RESULT=%ERRORLEVEL%"
echo.
if "%HOSC_INSTALL_RESULT%"=="0" (echo Completed successfully.) else (echo Installation did not complete. Read the error above.)
pause
exit /b %HOSC_INSTALL_RESULT%
