#!/bin/bash

# orange - G
# yellow - P
# blue   - C,A,V,L,I,M,F,W
# green  - S,T,N,Q
# red    - D,E,R,K
# cyan 	 - H,Y
# grey	 - B,X,U,Z,J,O

A3M=$1
GRM=$2

LEN=$(head -n 2 $A3M | tail -1 | awk '{print length($1)}')
ROW=$(grep -c "^>" $A3M)

alifile=$(mktemp /tmp/XXXXXX.ali)
grmfile=$(mktemp /tmp/XXXXXX.grm)
gplfile=$(mktemp /tmp/XXXXXX.gpl)
gapfile=$(mktemp /tmp/XXXXXX.gap)

grep -v "^>" $A3M | \
	sed "s/[a-z]//g; s/-/ 0/g; s/G/ 1/g; s/P/ 2/g" | \
	sed "s/[C,A,V,L,I,M,F,W]/ 3/g; s/[S,T,N,Q]/ 4/g; s/[D,E,R,K]/ 5/g" | \
	sed "s/[H,Y]/ 6/g; s/[B,X,U,Z,J,O]/ 7/g"  > $alifile

grep -v "^#" $GRM | awk '{if($3-$1>2) {print $1, $3, $5}}' | sort -k3 -g -r | \
	head -n $((LEN*3/2)) | sort -k3 -g > $grmfile

head -n 4 $GRM | tail -1 | \
	awk -v L="$LEN" '{
		split($2,c,"")
		for(i=1;i<=length($2);i++) {
			if(c[i]=="-") {
				print "set object rect from "i-1.5",-1 to "i-0.5","L" fc rgb \"#EEEEEE\" fillstyle solid 1.0 border";
				print "set object rect from -1,"i-1.5" to "L","i-0.5" fc rgb \"#EEEEEE\" fillstyle solid 1.0 border";
			}
		}
	}' > $gapfile

cat << EOF > $gplfile

set terminal pdf size 3.0,4.0 dashed dashlength 0.5
set output "contacts.pdf"

set multiplot

set tmargin at screen 0.90
set bmargin at screen 0.65
set rmargin at screen 0.85
set lmargin at screen 0.15

set yrange [$ROW:-1]
set xrange [-1:$LEN]
set cbrange [0:7]

set xtics out nomirror scale 0.5 font "Helvetica,10" format "%g" offset 0.0,0.5
set ytics out nomirror scale 0.5 font "Helvetica,10" format "%g" offset 0.6,0.0

set xlabel "" font "Helvetica,12" offset 0.0,1.0
set ylabel "" font "Helvetica,12" offset -1.0,0.0

set palette defined (0 "white", 1 "orange", 2 "yellow", 3 "blue", 4 "green", 5 "red", 6 "cyan", 7 "grey")
unset key
unset colorbox
set grid front

set obj 1 rect from screen 0.20,0.975 to screen 0.23,0.955 fc rgb "orange" fillstyle solid 1.0 noborder
set obj 2 rect from screen 0.20,0.940 to screen 0.23,0.920 fc rgb "yellow" fillstyle solid 1.0 noborder
set obj 3 rect from screen 0.30,0.975 to screen 0.33,0.955 fc rgb "blue" fillstyle solid 1.0 noborder
set obj 4 rect from screen 0.30,0.940 to screen 0.33,0.920 fc rgb "green" fillstyle solid 1.0 noborder
set obj 5 rect from screen 0.65,0.975 to screen 0.68,0.955 fc rgb "red" fillstyle solid 1.0 noborder
set obj 6 rect from screen 0.50,0.940 to screen 0.52,0.920 fc rgb "cyan" fillstyle solid 1.0 noborder
set obj 7 rect from screen 0.65,0.940 to screen 0.68,0.920 fc rgb "grey" fillstyle solid 1.0 noborder

set label 1 "G" at screen 0.24,0.966 font "Helvetica,8"
set label 2 "P" at screen 0.24,0.931 font "Helvetica,8"
set label 3 "C,A,V,L,I,M,F,W" at screen 0.34,0.966 font "Helvetica,8"
set label 4 "S,T,N,Q" at screen 0.34,0.931 font "Helvetica,8"
set label 5 "D,E,R,K" at screen 0.69,0.966 font "Helvetica,8"
set label 6 "H,Y" at screen 0.53,0.931 font "Helvetica,8"
set label 7 "unknown" at screen 0.69,0.931 font "Helvetica,8"

plot "$alifile" matrix w image

set tmargin at screen 0.59
set bmargin at screen 0.06

set yrange [$LEN:-1]
set xrange [-1:$LEN]
set cbrange [0:1]

set xtics out mirror

set style line 2604 lw .5 lc rgb "black"
set colorbox border 2604
set cbtics out nomirror scale 0.5 format "%.1f" font "Helvetica,8" offset -0.5,0.0

set palette defined (0 "white", 0.4 "white", 0.7 "blue", 1.0 "red")

load "$gapfile"
plot "$grmfile" u 1:2:3 palette pt 7 ps 0.25, \
	"" u 2:1:3 palette pt 7 ps 0.25, \
	x lw 1 lc rgb "grey"

EOF

gnuplot $gplfile

rm $gplfile
rm $alifile
rm $grmfile
rm $gapfile

