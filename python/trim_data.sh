#! /usr/bin/bash
INPUT=$1
EXTDMP=".dump"
EXTGRO=".gro"
EXTTGRO=".tgro"
EXTBOX=".box"
EXTTRIM=".trim"
GROOUT=$INPUT$EXTTGRO
BOXOUT=$INPUT$EXTBOX
TRIMOUT=$INPUT$EXTTRIM

sed -E 's/\s+/ /g' $INPUT$EXTGRO | sed -E 's/^\s//g' > $GROOUT
echo "Saved gro data to $GROOUT"

cat $INPUT$EXTDMP | grep -E " x\[" | grep -oP "{(.*)}" | sed "s/[{}]//g" | sed "s/\s//g" > $TRIMOUT
echo "Saved trim data to $TRIMOUT"

grep -E "box\[" $INPUT$EXTDMP | sed -E -e 's/^.*\{//' -e 's/0.0{5}e\+00//g' -e 's/[ ,}]//g' > $BOXOUT
echo "Saved box data to $BOXOUT"
