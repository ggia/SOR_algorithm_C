#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* C supporting bool types */
#include <stdbool.h>
/* Global variables */
double e = 0.0000005;	/* accurancy (error): 1/2 * 10^-6   */
int MAXITER = 5000;		/* max number of solving iterations */

double *coef; 			/* a, b, c */
double *lysh;
double *d;
int size;			/* Size of Matrix A */

int SOR(double, double, int);

int main(int argc, char **argv) {

	int a;		/* Παράμετρος α */
	int b;		/* Παράμετρος β */

	/* Αν έχουμε παράμετρους από γραμμή εντολών.. */
	if (argc == 4) {
		a = atoi(argv[1]);
		b = atoi(argv[2]);
		size = atoi(argv[3]);
	}

	else {
		printf("Παράμετροι από γραμμή εντολών: %s α β n \n\n", argv[0]);
		printf("Π.χ. τρέχουμε \"%s 1 2 100\" για το α=1, β=2 και n=100 \n\n", argv[0]);
		printf("Δώσε τιμές α, β και και της διάστασης A: ");
		scanf("%d%d%d", &a, &b, &size);
	}

	printf("\nΈγινε εισαγωγή a = %d, β = %d, n = %d \n\n", a, b, size);


	lysh = (double *)malloc(sizeof(double)*size);
	d = (double *)malloc(sizeof(double)*size);					/* πίνακας d */

	int minIter;

	coef = (double *)malloc(3 * sizeof(double));

	int i;
	for (i = 0; i < size; i++) {
		coef[1] = 4;		/* διαγώνιος 4 */
		coef[0] = -a;		/* a=-α		   */
		coef[2] = -b;		/* b=-β		   */
	}

	/* θεωρώ την λύση X=(1,1,..,1) και ορίζω τιμές στο d = Ax */
	d[0] = 4 - b;
	d[size - 1] = 4 - a;
	for (i = 1; i <= size - 2; i++)
		d[i] = 4 - a - b;

	double w, wb;

	minIter = MAXITER;
	for (w = 0.1; w <= 2; w += 0.1) {
		int iter = SOR(w, e, MAXITER);
		if (iter > 0 && iter < minIter) {
			minIter = iter;
			wb = w;
		}
		printf("Με w= %f έτρεξε %d φορές \n", w, iter);
	}

	/* Έχει βρεθεί λύση; */
	if (minIter != MAXITER) {
		/* Τρέχει η υπορουτίνα SOR */
		SOR(wb, e, MAXITER);
		/* Εκτύπωσε λύση */
		printf("\n\nΑκολουθούν ο πίνακας x και d από το σύστημα Ax=d\n");
		printf("\nx = (");
		for (i = 0; i < size - 1; i++)
			printf("%f, ", lysh[i]);
		printf("%f)\n", lysh[size - 1]);
		printf("\n\twb = %f, ελάχιστος αριθμός επαναλήψεων %d \n", wb, minIter);
		printf("\nd = (");
		for (i = 0; i < size - 1; i++)
			printf("%f, ", d[i]);
		printf("%f)\n", d[size - 1]);

	}
	else
		printf("Δεν βρέθηκαν λύσεις μετά από %d επαναλήψεις!\n", MAXITER);


	/* αποδέσμευση πόρων */
	free((void**)coef);
	free((void*)d);
	free((void*)lysh);

	return EXIT_SUCCESS;
}

/* Συνάρτηση που υλοποιεί την μέθοδο SOR */
int SOR(double w, double e, int maxIter) {

	int Epanalhpseis = 0;
	int i;

	double *xPalio = (double *)malloc(sizeof(double)*size);
	double *xTwra = (double *)malloc(sizeof(double)*size);
	double s1 = 0, s2 = 0, s3 = 0;

	bool solved = false; 	/* false */

	for (i = 0; i < size; i++)
		xPalio[i] = d[i];

	while (Epanalhpseis <= maxIter)
	{
		for (i = 0; i < size; i++)
		{
			s1 = s2 = s3 = 0;

			if (i > 0) {
				s1 = -coef[0] / coef[1] * xTwra[i - 1];
				s2 = -coef[0] / coef[1] * xPalio[i - 1];
			}

			if (i < size - 1) {
				s3 = -coef[2] / coef[1] * xPalio[i + 1];
			}

			/* SOR συνάρτηση */
			xTwra[i] = (1 - w) * xPalio[i] + w * (s1 + s3 + (double)d[i] / coef[1]);
		}

		Epanalhpseis++;

		bool syglish = true;
		/* έλεγχος επιθυμητής ακρίβειας προσέγγισης  */
		for (i = 0; i < size; i++) {
			double diff = xTwra[i] - xPalio[i];
			if (diff<0)
				diff = -diff;
			if (diff >= e)
				syglish = false;
		}

		if (syglish) {
			solved = true;
			break;
		}

		for (i = 0; i < size; i++)
			xPalio[i] = xTwra[i];
	}


	if (solved) {
		for (i = 0; i < size; i++)
		{
			if (!  (xTwra[i] == xTwra[i]) && 
					xTwra[i] <= DBL_MAX &&
					xTwra[i] >= -DBL_MAX) {
				solved = false;
				break;
			}
		}

		if (solved)
		{
			for (i = 0; i < size; i++) {
				lysh[i] = xTwra[i];
			}
		}
	}

	free((void*)xPalio);
	free((void*)xTwra);

	if (solved)
		return Epanalhpseis;
	else
		return -1;
}

