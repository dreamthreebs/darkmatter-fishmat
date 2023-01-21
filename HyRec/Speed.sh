rm tmp.dat
date > tmp.dat
for ((k=1;k<100;k=k+1))
do
./hyrec<input.dat>output.dat
done
date >> tmp.dat
