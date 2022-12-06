all: nbody_it generate_nbody

options_dir = ./opts/
nbody_it: 
	cc -o nbody_it ./src/nbody_it.c ./src/lib/*.c -lm
	./nbody_it -p ./opts/default -x

generate_nbody:
	cc -o generate_nbody ./src/generate_nbody.c ./src/lib/*.c -lm

$(options_dir):
	mkdir $@

plots: nbody_it
	./nbody_it -m leap -h 0.05 -e 0.5 -s 9000 -n 9 -O data/5_ecc_leap_out.csv
	./nbody_it -m rk4 -h 0.05 -e 0.5 -s 9000 -n 9 -O data/5_ecc_rk4_out.csv
	./nbody_it -m leap -h 0.003 -e 0.9 -s 150000 -n 150 -O data/9_ecc_leap_out.csv
	./nbody_it -m rk4 -h 0.003 -e 0.9 -s 150000 -n 150 -O data/9_ecc_rk4_out.csv
	python ./src/plots.py

clean:
	rm -f nbody_it
	rm -f generate_nbody

refresh: clean all

backup:
	rm -rf ../.term_proj_backup/*
	cp -r * ../.term_proj_backup
