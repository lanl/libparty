rsync -aP darwin:~/libparty/*.txt .
tail -n +2 nnaos.txt > nnaos.tmp
tail -n +2 nnaosoa.txt > nnaosoa.tmp
tail -n +2 nnsoa.txt > nnsoa.tmp
tail -n +2 aos.txt > aos.tmp

mv nnaos.tmp nnaos.txt
mv nnaosoa.tmp nnaosoa.txt
mv nnsoa.tmp nnsoa.txt
mv aos.tmp aos.txt

gnuplot nnompscaling.plot
gnuplot ompscaling.plot
gnuplot nnscaling.plot
gnuplot stepscaling.plot
convert nnomp.svg nnomp.png
convert stepomp.svg stepomp.png
convert nnn.svg nnn.png
convert stepn.svg stepn.png

display nnomp.png &
display stepomp.png &
display nnn.png &
display stepn.png &


