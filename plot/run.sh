#!/bin/sh

#aa=8
#bb=19
for aa in 0 1 2 3 4 5 6 7 8 9 10
do
for bb in 11 12 13 14 15 16 17 18 19
do
sed -i "23s/^.*.$/hgt_anom = dim_avg_n_Wrap(Hgt_anom($aa:$bb,:,:), 0)/g" wave-P1.ncl
sed -i "24s/^.*.$/hgt_clim = dim_avg_n_Wrap(Hgt_clim($aa:$bb,:,:), 0)/g" wave-P1.ncl
sed -i "167s/^.*.$/wks = gsn_open_wks(\"eps\", \"wave-\"+case+\"-$aa-$bb\")/g" wave-P1.ncl

sed -i "23s/^.*.$/hgt_anom = dim_avg_n_Wrap(Hgt_anom($aa:$bb,:,:), 0)/g" wave-P2.ncl
sed -i "24s/^.*.$/hgt_clim = dim_avg_n_Wrap(Hgt_clim($aa:$bb,:,:), 0)/g" wave-P2.ncl
sed -i "167s/^.*.$/wks = gsn_open_wks(\"eps\", \"wave-\"+case+\"-$aa-$bb\")/g" wave-P2.ncl

ncl wave-P1.ncl
ncl wave-P2.ncl
done
done

for ii in `ls *.eps`
do
	convert -density 800 -trim $ii ${ii/eps/png}
	echo $ii
done

rm -rf *.eps
