int trmIndex(int, int);

double ed(double*, int, int, int);

double stress(double, double, int, int);

int neighbours(double, int, double, int*, int*, int*);

int cmp_ivwrapper(const void*, const void*);

int nearest_neighbours(double*, int, int);

double CkNNStochasticProximityEmbedding(
        double*,
        double*,
        int,
        int,
        int,
        double,
        double,
        int,
        int,
        int);

double CkNeighboursStochasticProximityEmbedding(
	double*,
        double*,
        double,
        int,
        int,
        int,
        double,
        double,
        int,
        int);

double CStochasticProximityEmbedding(
        double*,
        double*,
        double,
        int,
        int,
        double,
        double,
        int,
        int,
        int);




