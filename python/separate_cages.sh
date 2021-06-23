#! /usr/bin/bash
INPUT=$1
EXTIN=".out"
EXTSMA="-small.dat"
EXTLAR="-large.dat"
EXTINT="-inter.dat"
EXTIRR="-irreg.dat"

grep -E "Small cages" $INPUT$EXTIN | sed -E 's/[a-z]//g' | sed 's/[A-Z]//g' | sed -E 's/\s\:\s//g' > $INPUT$EXTSMA
grep -E "Large cages" $INPUT$EXTIN | sed -E 's/[a-z]//g' | sed 's/[A-Z]//g' | sed -E 's/\s\:\s//g' > $INPUT$EXTLAR
grep -E "Interface cages" $INPUT$EXTIN | sed -E 's/[a-z]//g' | sed 's/[A-Z]//g' | sed -E 's/\s\:\s//g' > $INPUT$EXTINT
grep -E "Irregular cages" $INPUT$EXTIN | sed -E 's/[a-z]//g' | sed 's/[A-Z]//g' | sed -E 's/\s\:\s//g' > $INPUT$EXTIRR

echo "Done"
