#First run
#python format_ensemble.py 12 3
#python run_ensemble.py 1 1 1 0
#mpirun -n 4 python load_synthetics.py 1
#mpirun -n 4 python forward_parallel.py 1
#python run_ensemble.py 2 0 1 0
#mpirun -n 4 python load_synthetics.py 2
#python run_ensemble.py 2 0 0 1


#python format_ensemble.py Illapel 10 3
#python run_ensemble.py Illapel 1 1 1 0
#mpirun -n 12 python forward_parallel.py Illapel 1
#for (( i=2; i<=10; i++));
#do
#	echo $i 
#	python run_ensemble.py Illapel $i 0 1 0
#	mpirun -n 12 python forward_parallel.py Illapel $i
#	rm -rf Illapel/100_samples/$i/GFs/dynamic
#done

#python format_ensemble.py Tohoku 10 2
#python run_ensemble.py Tohoku 1 1 1 0
#mpirun -n 12 python forward_parallel.py Tohoku 1
#for (( i=2; i<=10; i++));
#do
#	echo $i 
#	python run_ensemble.py Tohoku $i 0 1 0
#	mpirun -n 12 python forward_parallel.py Tohoku $i
#	rm -rf Tohoku/100_samples/$i/GFs/dynamic
#done

#python format_ensemble.py Gorkha 10 2
#python run_ensemble.py Gorkha 1 1 1 0
#mpirun -n 12 python forward_parallel.py Gorkha 1
#for (( i=2; i<=10; i++));
#do
#	echo $i 
#	python run_ensemble.py Gorkha $i 0 1 0
#	mpirun -n 12 python forward_parallel.py Gorkha $i
#	rm -rf Gorkha/100_samples/$i/GFs/dynamic
#done

#python format_ensemble.py Pedernales 10 2
#python run_ensemble.py Pedernales 1 1 1 0
#mpirun -n 12 python forward_parallel.py Pedernales 1
#for (( i=2; i<=10; i++));
#do
#	echo $i 
#	python run_ensemble.py Pedernales $i 0 1 0
#	mpirun -n 12 python forward_parallel.py Pedernales $i
#	rm -rf Pedernales/100_samples/$i/GFs/dynamic
#done

name=Pedernales
python utils/format_ensemble.py $name 10 1
python utils/run_ensemble.py $name 1 1 1 0
mpirun -n 4 python utils/forward_parallel.py $name 1
for (( i=2; i<=10; i++));
do
	echo $i 
	python utils/run_ensemble.py $name $i 0 1 0
	mpirun -n 4 python utils/forward_parallel.py $name $i
	rm -rf Dynamic_Simulations/$name/100_samples/$i/GFs/dynamic
done

