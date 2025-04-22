load "paddim_v3_rotation.dat"
set terminal postscript enh col "Times-Roman,14"
set output "subplots.eps"
set multiplot layout 2,2 columnsfirst margins 0.1,0.99,0.02,0.99 spacing 0.12, .15

# first plot
set ylabel "|{/Bold u}_h|_{rms}" rotate by 0
set xlabel "1/Ro"
set log xy
set xrange [0.4:10]
set key top left
set key spacing 1.3 font ",12"
plot sample\
 [i=1:5] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:8] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1", \
 [i=13:15] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle, \
 [1:4] 2*x dt 2 lw 2 lc rgb "dark-blue" title "2Ro^{-1}"

# -------------------------------------------------------------

# second plot

set xlabel "1/Ro"
set ylabel "{/Symbol c}" rotate by 0
set format y "10^{%T}"
set log xy
set xrange [0.4:10]
#set yrange [0.01:1]
set key bottom left spacing 1.3 font ",12"

plot sample\
 [i=1:5] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:8] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1",\
 [i=13:15] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle



# -------------------------------------------------------------

# third plot

set xlabel "1/Ro"
set ylabel "w_{rms}" rotate by 0
set key  bottom left spacing 1.3 font ",12"
unset format y
set log xy
set xrange [0.4:10]
#set yrange [0.1:1]

plot sample \
 [i=1:5] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:8] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1",\
 [i=13:15] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle



# -------------------------------------------------------------

# fourth plot

set xlabel "1/Ro"
set ylabel "{/Symbol h}" rotate by 0
set key bottom left spacing 1.3 font ",12"
set log xy
set xrange [0.4:10]
#set yrange [0.1:1]

plot sample\
 [i=1:5] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:8] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1",\
 [i=13:15] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle


