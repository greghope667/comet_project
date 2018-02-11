#!/bin/bash
usage() {
echo "Usage:
$0 [-h] [-sk] [-t THREADS] [-o OF] -w WORKDIR file1.tgz [file2.tgz ...]

Runs analysis on light curves in target compressed archives.

Parameters:
    -h  Display this text.
    -k  Keeps temporary files on completion.
        (do not use with multiple archives)
    -o  File to output results to (default ./output.txt).
    -s  Split output to multiple files.
    -t  Number of threads to use (default is 1).
    -w  Working directory to store temporary extracted files."
exit 1
}

format_time() {
    local seconds=$1
    local h=$(( $seconds / 3600 ))
    local m=$(( $seconds % 3600 / 60 ))
    local s=$(( $seconds % 60 ))
    printf "%02d:%02d:%02d" $h $m $s
}

# Default options
threads=1
outputfile="output.txt"
workdir=""
keepfiles=false
splitout=false

while getopts ':hskt:o:w:' opt; do
    case $opt in
        h)
            usage
            ;;
        t)
            threads=${OPTARG}
            ;;
        o)
            outputfile=${OPTARG}
            ;;
        w)
            workdir=${OPTARG}
            ;;
        k)
            keepfiles=true
            ;;
        s)
            splitout=true
            ;;
        *)
            usage
            ;;
    esac
done

if [ ! -d "$workdir" ]; then
    echo "Working directory \'$workdir\' does not exist, exiting."
    exit
fi

shift $((OPTIND-1))

for file in "$@"; do
    echo "Current file: $file"
    echo "Checking disk space requirements"
    size=$(tar tzvf $file | awk '{s+=$3} END{print s}')
    space=$(df -P --block-size=1 $workdir | tail -n 1 | awk '{print $4}')

    if (($size > $space)); then
        echo "Not enough space to extract file, skipping"
        continue
    fi

    echo "Extracting archive"
    tar xzf $file -C $workdir

    file_count=$( find $workdir -name '*.fits' | wc -l )
    echo "$file_count files extracted"

    if [ -f $workdir/out.txt ]; then
        rm $workdir/out.txt
    fi

    ./batch_analyse.py -t $threads -o $workdir/out.txt $workdir &

    start_time=$SECONDS

    while [ ! -f $workdir/out.txt ]
    do
        sleep 1
    done

    while true; do
        count=$(wc -l $workdir/out.txt | awk '{print $1}')
        elapsed_time=$(( $SECONDS - $start_time ))
        remaining_count=$(( $file_count - $count ))
        remaining_time=$(( $remaining_count * $elapsed_time / $count ))

        printf "\rProcessed: $count   ETA: $( format_time $remaining_time )"

        if (($count == $file_count)); then
            break
        fi
        sleep 1
    done

    printf "\n"
    sleep 1

    if $splitout; then
        suffix="$(basename $file .tgz)"
    else
        suffix=""
    fi

    cat $workdir/out.txt >> $outputfile$suffix

    if ! $keepfiles; then
        echo "Cleaning temporary files"
        find $workdir -name '*.fits' | xargs rm
        rm $workdir/out.txt
    fi

    echo "Elapsed time: $( format_time $(( $SECONDS - $start_time )) )"
done
