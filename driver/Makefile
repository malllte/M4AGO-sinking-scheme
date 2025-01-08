default: run

reconfigure:
	meson --reconfigure ../build

setup:
	meson setup ../build

compile: setup
	meson compile -C ../build

run: compile
	meson test -t 60 -C ../build  "run M4AGO"

clean:
	rm -rf ../build/*
