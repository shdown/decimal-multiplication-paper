#!/usr/bin/env bash

# (c) 2020 shdown
# This code is licensed under MIT license (see LICENSE.MIT for details)

# This script expects the following binaries to be present:
#   * './bench_mpdec': basically just compiled 'bench_mpdec.c' (link it with '-lmpdec').
#   * './bench_fast': compiled from 'bench.c', 'fft.c', 'fft_utils.c', with
#                         #define FFT_USE_BUILTIN_UNPREDICTABLE() 1
#                     in 'fft.h'.
#   * './bench_slow': compiled from 'bench.c', 'fft.c', 'fft_utils.c', with
#                         #define FFT_USE_BUILTIN_UNPREDICTABLE() 0
#                     in 'fft.h'.
#
# For each N, it outputs three columns:
#   N  ratio_fast  ratio_slow
#
# Redirect output of this script to a file, then post-process that file with
#  awk '{print log($1)/log(2), $2, $3}' raw.txt > cooked.txt

set -e

bench() {
    local t_0f t_0s t_d
    t_0f=$(./bench_fast 0 "$1" "$2") || return $?
    t_0s=$(./bench_slow 0 "$1" "$2") || return $?
    t_d=$(./bench_mpdec "$1" "$2") || return $?

    ratio1=$(echo "scale=4; $t_d / $t_0f" | bc) || return $?
    ratio2=$(echo "scale=4; $t_d / $t_0s" | bc) || return $?

    echo "$1 $ratio1 $ratio2"
    echo >&2
    echo >&2 "~~~~~> ratios: $ratio1, $ratio2"
    echo >&2
}

get_niters_for_n() {
    echo $(( 80000000 / $1 ))
}

best_case_ours() {
    local x
    for (( x = 128; x <= 4096; x *= 2 )); do
        echo "$(( 17 * x ))"
        echo "$(( 17 * (x/2*3) ))"
    done

    echo "$(( 17 * 8192 ))"

    echo 262144
    echo 393216

    echo 524288
    echo 786432

    echo 1048576
    echo 1572864

    echo 2097152
    echo 3145728

    echo 4194304
    echo 6291456

    echo 8388608
    echo 12582912

    echo 15728640
    echo 23592960
}

best_case_mpdec() {
    local x
    for (( x = 128; x <= 1048576; x *= 2 )); do
        echo "$(( 19 * (x/2*3) ))"
    done
}

weird_case_mpdec() {
    local x
    for (( x = 128; x <= 1048576; x *= 2 )); do
        echo "$(( 19 * x ))"
    done
}

ns=$({ best_case_ours; best_case_mpdec; weird_case_mpdec; } | sort -n)
prev=0
for n in $ns; do
    if (( prev != 0 )); then
        if (( n - prev < 2 )); then
            echo >&2 "E: staircase collapse."
            exit 1
        fi
        avg=$(( (n + prev) / 2 ))
        bench "$avg" "$(get_niters_for_n "$avg")"
    fi
    prev=$n
done

last=$(( prev + 1 ))
bench "$last" "$(get_niters_for_n "$last")"
