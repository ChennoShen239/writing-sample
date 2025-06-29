@echo off
echo Starting compilation of Sovereign Default Model...

REM Setup Intel oneAPI environment variables
echo Setting up Intel oneAPI environment...
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 > nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo Warning: Could not set Intel oneAPI environment variables, this may affect program execution
)

REM Create bin directory if it doesn't exist
if not exist bin mkdir bin
echo Created bin directory (if it didn't exist before)

REM Compile source files
echo Compiling sim.f90...
ifx /O3 /Qopenmp src\sim.f90 /c

echo Compiling NL.f90...
ifx /O3 /Qopenmp src\NL.f90 /c

echo Compiling defMod.f90...
ifx /O3 /Qopenmp src\defMod.f90 /c

echo Compiling defaultModel.f90 and linking all object files...
ifx /O3 /Qopenmp /Qmkl sim.obj NL.obj defMod.obj src\defaultModel.f90 /Fe:bin\defaultModel.exe

if %ERRORLEVEL% equ 0 (
    echo Compilation successful! Executable file generated: bin\defaultModel.exe
    
    REM Copy runtime libraries to bin directory to ensure the program can run
    echo Copying necessary runtime libraries...
    if exist "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler\libiomp5md.dll" (
        copy "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler\libiomp5md.dll" bin\ > nul
        echo Copied libiomp5md.dll to bin directory
    ) else (
        echo Warning: Could not find libiomp5md.dll, the program may not run
    )
    
    echo.
    echo To run the program, use the command: bin\defaultModel.exe
    echo After running, you can use MATLAB to analyze the results in the results directory
) else (
    echo Compilation failed, please check the error messages
)

echo.
echo Compilation process completed 