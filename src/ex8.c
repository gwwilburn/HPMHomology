#include <stdio.h>

int main(int argc, char *argv[])
{
	int i = 0;

	if (argc == 2) {
		printf("You only have one argument. :( \n");
	} else if (argc > 2 && argc < 5) {
		printf("Here are your arguments:\n");

		for (i=0; i < argc; i++) {
			printf("%d: %s \n", i, argv[i]);
		}
	} else {
		printf("You have too many arguments.\n\n");

	}


	return 0;

}


