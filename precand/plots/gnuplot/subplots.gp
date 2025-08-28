load "paddim_v3_rotation.dat"
load "highres_rotation.dat"
load "newIC_rotation.dat"
set terminal postscript enh col "Times-Roman,14" size 8in,12in
set output "highres_subplots.eps"
set multiplot layout 3,2 columnsfirst #margins 0.1,0.99,0.02,0.99 spacing 0.12, .15

# first plot
set ylabel "|{/Bold u}_h|_{rms}" rotate by 0
set xlabel "1/Ro"
set log xy
set xrange [0.4:10]
set key top left
set key bottom right spacing 1.3 font ",8"
plot sample\
 [i=1:1] '+' using (invRo_hr[i]):(sqrt(urms_hr[i]**2 + vrms_hr[i]**2)):(sqrt(uerr_hr[i]**2 + verr_hr[i]**2)) w yerrorbars pt 9 ps 2 lw 2 lc rgb "blue" title "TG Re=1000, Fr=0.18",\
 [i=2:2] '+' using (invRo_hr[i]):(sqrt(urms_hr[i]**2 + vrms_hr[i]**2)):(sqrt(uerr_hr[i]**2 + verr_hr[i]**2)) w yerrorbars pt 9 ps 2 lw 2 lc rgb "red" title "TG Re=1000, Fr=0.1 ",\
 [i=1:2] '+' using (invRo_ic[i]):(sqrt(urms_ic[i]**2 + vrms_ic[i]**2)):(sqrt(uerr_ic[i]**2 + verr_ic[i]**2)) w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" title "TG Re= 600, Fr=0.18",\
 [i=3:4] '+' using (invRo_ic[i]):(sqrt(urms_ic[i]**2 + vrms_ic[i]**2)):(sqrt(uerr_ic[i]**2 + verr_ic[i]**2)) w yerrorbars pt 5 ps 2 lw 2 lc rgb "red" title "TG Re= 600,  Fr=0.1 ",\
 [i=1:5] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Re = 600, Fr=0.18",\
 [i=6:8] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Re= 600, Fr=0.1 ", \
 [i=13:15] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle, \
 [1:4] 2*x dt 2 lw 2 lc rgb "dark-blue" title "2Ro^{-1}"

# -------------------------------------------------------------

# second plot

unset yrange
set xlabel "1/Ro"
set ylabel "{/Symbol c}" rotate by 0
set format y "10^{%T}"
set log xy
set xrange [0.4:10]
#set yrange [0.01:1]
set key bottom left spacing 1.3 font ",8"

plot sample\
 [i=1:1] '+' using (invRo_hr[i]):(Dtemp_hr[i]):(dterr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "blue" title "TG Re=1000, Fr=0.18",\
 [i=2:2] '+' using (invRo_hr[i]):(Dtemp_hr[i]):(dterr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "red" title "TG Re=1000, Fr=0.1 ",\
 [i=1:2] '+' using (invRo_ic[i]):(Dtemp_ic[i]):(dterr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" title "TG Re= 600, Fr=0.18",\
 [i=3:4] '+' using (invRo_ic[i]):(Dtemp_ic[i]):(dterr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "red" title "TG Re= 600, Fr=0.1 ",\
 [i=1:5] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "   Re= 600, Fr=0.18",\
 [i=6:8] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "   Re= 600, Fr=0.1 ",\
 [i=13:15] '+' using (invRo[i]):(Dtemp[i]):(dterr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle



# -------------------------------------------------------------

# third plot

set xlabel "1/Ro"
set ylabel "w_{rms}" rotate by 0
set key  bottom left spacing 1.3 font ",8"
unset format y
set log xy
set xrange [0.4:10]
#set yrange [0.1:1]

plot sample \
 [i=1:1] '+' using (invRo_hr[i]):(wrms_hr[i]):(werr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "blue" title "TG Re=1000, Fr=0.18",\
 [i=2:2] '+' using (invRo_hr[i]):(wrms_hr[i]):(werr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "red" title "TG Re=1000, Fr=0.1 ",\
 [i=1:2] '+' using (invRo_ic[i]):(wrms_ic[i]):(werr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" title "TG Re= 600, Fr=0.18",\
 [i=3:4] '+' using (invRo_ic[i]):(wrms_ic[i]):(werr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "red" title "TG Re= 600, Fr=0.1 ",\
 [i=1:5] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "   Re= 600, Fr=0.18",\
 [i=6:8] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "   Re= 600, Fr=0.1 ",\
 [i=13:15] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle



# -------------------------------------------------------------

# fourth plot

set xlabel "1/Ro"
set ylabel "{/Symbol h}" rotate by 0
set key bottom left spacing 1.3 font ",8"
set log xy
set xrange [0.4:10]
#set yrange [0.1:1]

plot sample\
 [i=1:1] '+' using (invRo_hr[i]):(mix_hr[i]):(mixerr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "blue" title "TG Re=1000, Fr=0.18",\
 [i=2:2] '+' using (invRo_hr[i]):(mix_hr[i]):(mixerr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "red" title "TG Re=1000, Fr=0.1 ",\
 [i=1:2] '+' using (invRo_ic[i]):(mix_ic[i]):(mixerr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" title "TG Re= 600, Fr=0.18",\
 [i=3:4] '+' using (invRo_ic[i]):(mix_ic[i]):(mixerr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "red" title "TG Re= 600, Fr=0.1 ",\
 [i=1:5] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "   Re= 600, Fr=0.18",\
 [i=6:8] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "   Re= 600, Fr=0.1 ",\
 [i=13:15] '+' using (invRo[i]):(mix[i]):(mixerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle


# -------------------------------------------------------------

# fifth plot

set xlabel "1/Ro"
set ylabel "T" rotate by 0
set key bottom left spacing 1.3 font ",8"
set log xy
set xrange [0.4:10]
set yrange [0.01:0.2]

plot sample\
 [i=1:1] '+' using (invRo_hr[i]):(trms_hr[i]):(terr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "blue" title "TG Re=1000, Fr=0.18",\
 [i=2:2] '+' using (invRo_hr[i]):(trms_hr[i]):(terr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "red" title "TG Re=1000, Fr=0.1 ",\
 [i=1:2] '+' using (invRo_ic[i]):(trms_ic[i]):(terr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" title "TG Re= 600, Fr=0.18",\
 [i=3:4] '+' using (invRo_ic[i]):(trms_ic[i]):(terr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "red" title "TG Re= 600, Fr=0.1 ",\
 [i=1:5] '+' using (invRo[i]):(trms[i]):(terr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "   Re= 600, Fr=0.18",\
 [i=6:8] '+' using (invRo[i]):(trms[i]):(terr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(trms[i]):(terr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "   Re= 600, Fr=0.1 ",\
 [i=13:15] '+' using (invRo[i]):(trms[i]):(terr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle


# -------------------------------------------------------------

# sixth plot

set xlabel "1/Ro"
set ylabel "{/Symbol w}_{z}" rotate by 0
set key bottom left spacing 1.3 font ",8"
set log xy
set xrange [0.4:10]
set yrange [1:10]

plot sample\
 [i=1:1] '+' using (invRo_hr[i]):(wzrms_hr[i]):(wzerr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "blue" title "TG Re=1000, Fr=0.18",\
 [i=2:2] '+' using (invRo_hr[i]):(wzrms_hr[i]):(wzerr_hr[i]) w yerrorbars pt 9 ps 2 lw 2 lc rgb "red" title "TG Re=1000, Fr=0.1 ",\
 [i=1:2] '+' using (invRo_ic[i]):(wzrms_ic[i]):(wzerr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" title "TG Re= 600, Fr=0.18",\
 [i=3:4] '+' using (invRo_ic[i]):(wzrms_ic[i]):(wzerr_ic[i]) w yerrorbars pt 5 ps 2 lw 2 lc rgb "red" title "TG Re= 600, Fr=0.1 ",\
 [i=1:5] '+' using (invRo[i]):(wzrms[i]):(wzerr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "   Re= 600, Fr=0.18",\
 [i=6:8] '+' using (invRo[i]):(wzrms[i]):(wzerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(wzrms[i]):(wzerr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "   Re= 600, Fr=0.1 ",\
 [i=13:15] '+' using (invRo[i]):(wzrms[i]):(wzerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle


