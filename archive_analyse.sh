#!/bin/bash

usage() {
echo "Usage:
$0 [-h] [-adks] [-t THREADS] [-o OF] -w WORKDIR file1 [file2 ...]

Runs analysis on light curves in target compressed archives.

Parameters:
    -a  Split output for multiple archives, appending archive name to output.
    -d  Process directories of fits files, instead of archive files.
    -h  Display this text.
    -k  Keeps temporary files on completion.
        (do not use with multiple archives)
    -o  File to output results to (default ./output.txt).
    -q  Keep only points with SAP_QUALITY=0.
    -s  Skip disk space check for extraction.
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
checkspace=true
processarchives=true
q_flag=""

while getopts ':hadksqt:o:w:' opt; do
    case $opt in
        h)
            usage
            ;;
        t)
            [ "${OPTARG}" -gt 0 ] 2>/dev/null && threads=${OPTARG}
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
        a)
            splitout=true
            ;;
        s)
            checkspace=false
            ;;
        d)
            processarchives=false
            keepfiles=true
            ;;
        q)
            q_flag='-q'
            ;;
        *)
            usage
            ;;
    esac
done

if ( $processarchives && [ ! -d "$workdir" ]); then
    echo "Working directory '$workdir' does not exist, exiting."
    exit
fi

shift $((OPTIND-1))

for file in "$@"; do

    if $processarchives; then
        if [ -f $file ]; then
            echo "Current file: $file"
        else
            echo "Archive $file not found, skipping."
            continue
        fi

        if $checkspace; then
            echo "Checking disk space requirements."
            size=$(tar tzvf $file | awk '{s+=$3} END{print s}')
            space=$(df -P --block-size=1 $workdir |tail -n 1| awk '{print $4}')

            if (($size > $space)); then
                echo "Not enough space to extract file, skipping."
                continue
            fi
        fi

        echo "Extracting archive..."
        tar xzf $file -C $workdir
    else
        if [ -d $file ]; then
            echo "Current directory: $file"
            workdir=$file
        else
            echo "Directory $file not found, skipping."
            continue
        fi
    fi

    file_count=$( find $workdir -maxdepth 1 -name '*.fits' | wc -l )
    if $processarchives; then
        echo "$file_count files extracted."
    else
        echo "$file_count files found."
    fi

    if [ -f $workdir/out.txt ]; then
        rm $workdir/out.txt
    fi

    ./batch_analyse.py $q_flag -t $threads -o $workdir/out.txt $workdir &

    start_time=$SECONDS

    while [ ! -f $workdir/out.txt ]
    do
        sleep 1
    done

    count_prev=0
    time_since_prev=0

    while true; do
        count=$(wc -l $workdir/out.txt | awk '{print $1}')
        elapsed_time=$(( $SECONDS - $start_time ))
        remaining_count=$(( $file_count - $count ))
        remaining_time=$(( $remaining_count * $elapsed_time / $count ))

        printf "\rProcessed: $count   ETA: $( format_time $remaining_time )"

        if (($count == $file_count)); then
            printf "\n"
            break
        elif (($count > $count_prev)); then
            count_prev=$count
            time_since_prev=0
        elif (($time_since_prev > 60)); then
            printf "\nProcessing timed out, assumed finished.\n"
            break
        else
            time_since_prev=$(($time_since_prev + 1))
        fi
 
        sleep 1
    done

    sleep 1

    if $splitout; then
        suffix="$(basename $file .tgz)"
    else
        suffix=""
    fi

    cat $workdir/out.txt >> $outputfile$suffix

    if ! $keepfiles; then
        echo "Cleaning temporary files."
        find $workdir -maxdepth 1 -name '*.fits' -type f -delete
        rm $workdir/out.txt
    fi

    echo "Elapsed time: $( format_time $(( $SECONDS - $start_time )) )"
done

