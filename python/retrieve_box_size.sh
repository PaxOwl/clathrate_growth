#! /usr/bin/bash
INPUT=$1
EXTIN="dump"
EXTOUT="box"
OUTPUT=$(echo "data/$INPUT" | sed "s/$EXTIN/$EXTOUT/")

grep -E "box\[" $INPUT | sed -E -e 's/^.*\{//' -e 's/0.0{5}e\+00//g' -e 's/[ ,}]//g' > data/$OUTPUT
echo "Saved box data to $OUTPUT"