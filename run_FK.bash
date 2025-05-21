name=Pedernales
python utils/format_ensemble.py $name 10 1
python utils/run_ensemble.py $name 1 1 1 0
mpirun -n 10 python utils/forward_parallel.py $name 1
for (( i=2; i<=10; i++));
do
	echo $i 
	python utils/run_ensemble.py $name $i 0 1 0
	mpirun -n 10 python utils/forward_parallel.py $name $i
	rm -rf Dynamic_Simulations/$name/100_samples/$i/GFs/dynamic
done

