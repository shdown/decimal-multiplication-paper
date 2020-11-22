#!/usr/bin/env bash

# (c) 2020 shdown
# This code is licensed under MIT license (see LICENSE.MIT for details)

# This script expects './main' binary to be present.
# It should be compiled from the following files: 'main.c', 'fft.c', 'fft_utils.c'.

set -e

if (( $# != 1 )); then
    cat >&2 <<__EOF__
USAGE: $0 METHOD
  where METHOD must be either of: 1 4 6.
METHOD of 1 means test the "straigh" FFT (function 'fft()').
METHOD of 4 means test the four-step FFT (function 'fft_fourstep()').
METHOD of 6 means test the six-step FFT (function 'fft_sixstep()').
__EOF__
    exit 2
fi
METHOD=$1

trim_to_n() {
    local n=$1
    local q=$((n/1024))
    local r=$((n%1024))
    if (( q != 0 )); then
        dd status=none bs=1024 count=$q \
            || return $?
    fi
    if (( r != 0 )); then
        dd status=none bs=$r count=1 \
            || return $?
    fi
}

gen_n_digits() {
    local n=$1
    tr -cd '1-9' < /dev/urandom | trim_to_n 1 \
        || return $?
    tr -cd '0-9' < /dev/urandom | trim_to_n $(( n - 1 )) \
        || return $?
    echo
}

infile=/tmp/sudm_in.txt
outfile1=/tmp/sudm_out1.txt
outfile2=/tmp/sudm_out2.txt

test_for_n() {
    { gen_n_digits "$n"; gen_n_digits "$n"; } > "$infile" \
        || return $?
    ./main "$METHOD" < "$infile" > "$outfile1" \
        || return $?
    python3 -c 'x = int(input()); y = int(input()); print(x*y)' < "$infile" > "$outfile2" \
        || return $?
    cmp -l -- "$outfile1" "$outfile2" \
        || return $?
}

echo >&2 "~~~~~~~~> n < 512..."
for (( n = 1; n < 512; ++n )); do
    echo -n >&2 " $n"
    test_for_n $n
done

echo >&2 " OK"

for (( n = 512; n <= 262144; n *= 2 )); do
    echo >&2 "~~~~~~~~> n = $n"
    test_for_n $n
done

rm -f -- "$infile" "$outfile1" "$outfile2"
