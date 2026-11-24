@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\install_x64.ps1" %*
exit /b %ERRORLEVEL%
