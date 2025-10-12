@echo off
echo 开始编译主权债务违约模型...

REM 设置Intel oneAPI环境变量
echo 设置Intel oneAPI环境变量...
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 > nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo 警告: 无法设置Intel oneAPI环境变量，可能会影响程序运行
)

REM 创建bin目录（如果不存在）
if not exist bin mkdir bin
echo 已创建bin目录（如果之前不存在）

REM 编译各个源文件
echo 编译sim.f90...
ifx /O3 /Qopenmp src\sim.f90 /c

echo 编译NL.f90...
ifx /O3 /Qopenmp src\NL.f90 /c

echo 编译defMod.f90...
ifx /O3 /Qopenmp src\defMod.f90 /c

echo 编译defaultModel.f90并链接所有目标文件...
ifx /O3 /Qopenmp /Qmkl sim.obj NL.obj defMod.obj src\defaultModel.f90 /Fe:bin\defaultModel.exe

if %ERRORLEVEL% equ 0 (
    echo 编译成功！可执行文件已生成: bin\defaultModel.exe
    
    REM 复制运行时库到bin目录以确保程序可以运行
    echo 复制必要的运行时库...
    if exist "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler\libiomp5md.dll" (
        copy "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler\libiomp5md.dll" bin\ > nul
        echo 已复制libiomp5md.dll到bin目录
    ) else (
        echo 警告: 无法找到libiomp5md.dll，程序可能无法运行
    )
    
    echo.
    echo 要运行程序，请使用命令: bin\defaultModel.exe
    echo 运行完成后，可以使用MATLAB分析results目录中的结果
) else (
    echo 编译失败，请检查错误信息
)

echo.
echo 编译过程完成 