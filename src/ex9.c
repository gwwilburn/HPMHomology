#include <stdio.h>

int main(int argc, char *argv[])
{

	int i = 0;

	while (i < 25) {
		printf(" %d", i);
		i++;
	}

	//need this for the final newline
	printf("\n");

	while (i > 0) {
		printf(" %d", i);
		i--;
	}

	printf("\n");


	return 0;
}

