all: nbody_it generate_nbody

options_dir = ./opts/
nbody_it: 
	cc -o nbody_it ./src/nbody_it.c ./src/lib/*.c -lm
	./nbody_it -p ./opts/default -x

generate_nbody:
	cc -o generate_nbody ./src/generate_nbody.c ./src/lib/*.c -lm

$(options_dir):
	mkdir $@

clean:
	rm -f nbody_it
	rm -f generate_nbody

refresh: clean all

backup:
	rm -rf ../.term_proj_backup/*
	cp -r * ../.term_proj_backup
