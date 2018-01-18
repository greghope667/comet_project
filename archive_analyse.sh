#!/bin/bash
usage()
{
echo "(basename "$0") [-h] [-t THREADS] [-o OF] -w WORKDIR file1 [file2 ...]

Runs analysis on light curves in target zip files."
exit 1
}

# Default options
threads=1
outputfile="output.txt"
workdir=""

while getopts 'ht:o:w:' opt; do
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

    while [ ! -f $workdir/out.txt ]
    do
        sleep 1
    done

    while true; do
        count=$(wc -l $workdir/out.txt | awk '{print $1}')
        printf "\rProcessed: $count"
        if (($count == $file_count)); then
            break
        fi
        sleep 1
    done

    printf "\n"
    sleep 1

    cat $workdir/out.txt >> $outputfile

    echo "Cleaning temporary files"
    find $workdir -name '*.fits' | xargs rm
    rm $workdir/out.txt
done
