#!/bin/bash

# Change directory to the script's location to ensure relative paths work correctly.
cd "$(dirname "$0")" || exit

echo "Cleaning previous builds..."
rm -rf bin/* ./*.mod ./*.o

echo "Creating result directories..."
mkdir -p results/baseline
mkdir -p results/thetad100
mkdir -p results/theta10

echo "Building the model..."
# 编译模型（只需编译一次）
gfortran -Ofast -mcpu=apple-m1 -pipe -fopenmp -flto -Wall -pedantic -g -c src/NL.f90
gfortran -Ofast -mcpu=apple-m1 -pipe -fopenmp -flto -Wall -pedantic -g -c src/sim.f90
gfortran -Ofast -mcpu=apple-m1 -pipe -fopenmp -flto -Wall -pedantic -g -c src/defMod.f90
gfortran -Ofast -mcpu=apple-m1 -pipe -fopenmp -flto -Wall -pedantic -g -o bin/defaultModel *.o src/defaultModel.f90 -llapack -lblas -lm

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

# 函数：清理文件名中的空格
cleanup_filenames() {
    local dir=$1
    local case_name=$2
    
    echo "Cleaning and renaming files for $case_name..."
    
    # 进入目标目录
    cd "$dir"
    
    # 移除文件名中的前导空格并重命名
    for file in *; do
        # 移除前导空格
        new_name=$(echo "$file" | sed 's/^[[:space:]]*//')
        if [[ "$file" != "$new_name" ]]; then
            mv -- "$file" "$new_name"
            echo "  Trimmed: '$file' -> '$new_name'"
            file=$new_name # use new name for next step
        fi
    done

    # 专门重命名文件以匹配Matlab脚本的期望
    if [ -f "parameters.tab" ]; then
        mv parameters.tab par.dat
        echo "  Renamed: 'parameters.tab' -> 'par.dat'"
    fi
    if [ -f "sim.tab" ]; then
        mv sim.tab sim.dat
        echo "  Renamed: 'sim.tab' -> 'sim.dat'"
    fi

    # 返回原始目录
    cd - > /dev/null
}

echo "Running baseline case (thetaD=1)..."
export THETAD=1.0
export OUTDIR="./results/baseline/"
./bin/defaultModel

if [ $? -ne 0 ]; then
    echo "Baseline run failed!"
    exit 1
fi

# 清理baseline目录的文件名
cleanup_filenames "results/baseline" "baseline"

echo "Running high thetaD case (thetaD=100)..."
export THETAD=100.0
export OUTDIR="./results/thetad100/"
./bin/defaultModel

if [ $? -ne 0 ]; then
    echo "High thetaD run failed!"
    exit 1
fi

# 清理thetad100目录的文件名  
cleanup_filenames "results/thetad100" "thetad100"

echo "Running medium thetaD case (thetaD=50)..."
export THETAD=50.0
export OUTDIR="./results/theta10/"
./bin/defaultModel

if [ $? -ne 0 ]; then
    echo "Medium thetaD run failed!"
    exit 1
fi

# 清理theta10目录的文件名
cleanup_filenames "results/theta10" "theta10"

echo ""
echo "=== SIMULATIONS COMPLETED SUCCESSFULLY ==="
echo "✅ All three simulations completed"
echo "✅ File names cleaned (spaces removed)"
echo "✅ Ready for visualization"
echo ""
echo "Results saved to:"
echo "  - Baseline (thetaD=1): results/baseline/"
echo "  - High thetaD (thetaD=75): results/thetad100/"
echo "  - Medium thetaD (thetaD=10): results/theta10/"
echo ""
echo "🚀 Now running Matlab visualization..."

# 切换到results目录并运行Matlab
cd results
matlab -nodisplay -nosplash -nodesktop -r "try; main; fprintf('\n✅ All plots generated successfully!\n'); catch ME; fprintf('\n❌ Matlab visualization failed: %s\n', ME.message); end; exit;"

echo ""
echo "🎉 COMPLETE WORKFLOW FINISHED!"
echo "Check the results/ directory for output files and plots." 