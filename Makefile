detector: stat_detector.c main.c
	gcc -Wall -std=gnu11 -I. -c -g stat_detector.c
	gcc -Wall -std=gnu11 -I. -c -g main.c
	gcc -Wall -std=gnu11 -I. main.o stat_detector.o -lm -o detector
