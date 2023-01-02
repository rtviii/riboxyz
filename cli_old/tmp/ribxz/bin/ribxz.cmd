@echo off
setlocal enableextensions

if not "%RIBXZ_REDIRECTED%"=="1" if exist "%LOCALAPPDATA%\ribxz\client\bin\ribxz.cmd" (
  set RIBXZ_REDIRECTED=1
  "%LOCALAPPDATA%\ribxz\client\bin\ribxz.cmd" %*
  goto:EOF
)

if not defined RIBXZ_BINPATH set RIBXZ_BINPATH="%~dp0ribxz.cmd"
if exist "%~dp0..\bin\node.exe" (
  "%~dp0..\bin\node.exe" "%~dp0..\bin\run" %*
) else if exist "%LOCALAPPDATA%\oclif\node\node-14.17.6.exe" (
  "%LOCALAPPDATA%\oclif\node\node-14.17.6.exe" "%~dp0..\bin\run" %*
) else (
  node "%~dp0..\bin\run" %*
)
