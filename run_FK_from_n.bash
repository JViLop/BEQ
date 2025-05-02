#First run
#python format_ensemble.py 12 3
#python run_ensemble.py 1 1 1 0
#mpirun -n 4 python load_synthetics.py 1
#mpirun -n 4 python forward_parallel.py 1
#python run_ensemble.py 2 0 1 0
#mpirun -n 4 python load_synthetics.py 2
#python run_ensemble.py 2 0 0 1


python format_ensemble.py Illapel 100 3
#python run_ensemble.py Illapel 1 1 1 0
#mpirun -n 12 python forward_parallel.py Illapel 1
for (( i=61; i<=100; i++));
do
	echo $i 
	python run_ensemble.py Illapel $i 0 1 0
	mpirun -n 12 python forward_parallel.py Illapel $i
	rm -rf Illapel/100_samples/$i/GFs/dynamic
done

python format_ensemble.py Tohoku 100 2
#python run_ensemble.py Tohoku 1 1 1 0
#mpirun -n 12 python forward_parallel.py Tohoku 1
for (( i=61; i<=100; i++));
do
	echo $i 
	python run_ensemble.py Tohoku $i 0 1 0
	mpirun -n 12 python forward_parallel.py Tohoku $i
	rm -rf Tohoku/100_samples/$i/GFs/dynamic
done

python format_ensemble.py Gorkha 100 2
#python run_ensemble.py Gorkha 1 1 1 0
#mpirun -n 12 python forward_parallel.py Gorkha 1
for (( i=61; i<=100; i++));
do
	echo $i 
	python run_ensemble.py Gorkha $i 0 1 0
	mpirun -n 12 python forward_parallel.py Gorkha $i
	rm -rf Gorkha/100_samples/$i/GFs/dynamic
done

python format_ensemble.py Pedernales 100 2
#python run_ensemble.py Pedernales 1 1 1 0
#mpirun -n 12 python forward_parallel.py Pedernales 1
for (( i=61; i<=100; i++));
do
	echo $i 
	python run_ensemble.py Pedernales $i 0 1 0
	mpirun -n 12 python forward_parallel.py Pedernales $i
	rm -rf Pedernales/100_samples/$i/GFs/dynamic
done

python format_ensemble.py Iquique 100 3
#python run_ensemble.py Iquique 1 1 1 0
#mpirun -n 12 python forward_parallel.py Iquique 1
for (( i=61; i<=100; i++));
do
	echo $i 
	python run_ensemble.py Iquique $i 0 1 0
	mpirun -n 12 python forward_parallel.py Iquique $i
	rm -rf Iquique/100_samples/$i/GFs/dynamic
done

