#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <getopt.h>

using namespace std;


//Define global parameters
const int size = 100;				// Length of the square lattice
const int N = size * size;			// Number of spins
const double pi = M_PI;				// pi
const double thr_acceptance_rate = 0.5;		// threshold acceptance rate
int rs = 0;					// random seed

double r0 = 0.0;				// Temperature
double x = 1.0;					// detuning from degeneracy
double eps = 0.0;				// homogeneous strain
double u = 1.0;					// quartic interaction
double g = 10.0;				// quartic interaction
double lambda = 0.0;				// quartic interaction
double kappa = 0.0;				// quartic interaction
double lambda0 = 0.0;				// interaction to dislocation strain
double lambdax = 0.0;				// interaction to dislocation strain
double lambdaz = 0.0;				// interaction to dislocation strain
double K0 = 0.0;				// gradient term

int passes = 1000;				// Default # of passes
int pass = 0;					// pass number

double E = 0.0;					// Energy per site

double R1_state[N];				// Real part of Delta_1
double I1_state[N];				// Imaginary part of Delta_1
double R2_state[N];				// Real part of Delta_2
double I2_state[N];				// Imaginary part of Delta_2
double f[N];                                    // strain field


int nnup[N];					// nearest neighbor of site i above
int nndown[N];					// nearest neighbor of site i below
int nnleft[N];					// nearest neighbor of site i left
int nnright[N];					// nearest neighbor of site i right
int sweeporder[N];				// order in which sites are visited

double opening = 0.5;				// magnitude of Metropolis-Hastings fluctuations
double acceptance_rate = 1.0;			// acceptance rate of the finished pass
double new_acceptance_rate = 1.0;		// new acceptance rate of the running pass
double av_acceptance_rate = 0.0;		// average acceptance rate over all subpasses of one superpass

string printLastConf_Q = "False";		// flag defining whether the last configuration shall be printed

//Parses the input command line to redefine simulation parameters
void parse_input(int argc, char** argv)
{
	int c;
	while (1)
    {
    	static struct option long_options[] =
        {
        	/* These options don’t set a flag. We distinguish them by their indices. */
                {"r0",   required_argument, NULL, 't'},
                {"x",  required_argument, NULL, 'x'},
                {"eps",  required_argument, NULL, 'e'},
                {"u",  required_argument, NULL, 'y'},
                {"g",  required_argument, NULL, 'z'},
                {"lambda",  required_argument, NULL, 'j'},
                {"kappa",  required_argument, NULL, 'b'},
                {"lambda0",  required_argument, NULL, 'k'},
                {"lambdax",  required_argument, NULL, 'l'},
                {"lambdaz",  required_argument, NULL, 'm'},
                {"K0",  required_argument, NULL, 'i'},
        	{"passes",  required_argument, NULL, 'n'},
                {"printConfiguration",  required_argument, NULL, 's'},
        	{"rs",  required_argument, NULL, 'r'},
        	{NULL, 0, NULL, 0}
        };
    	/* getopt_long stores the option index here. */
    	int option_index = 0;

		c = getopt_long_only (argc, argv, "a:b:c:d:", long_options, NULL);

		/* Detect the end of the options. */
		if (c == -1)
        	break;
        	
        string soptarg = optarg;
		switch (c)
        {
        	case 0:
                        cout << "error" << endl;	break;
                case 't':
                        r0 = stod(optarg);		break;
                case 'x':
                        x = stod(optarg);               break;
                case 'e':
                        eps = stod(optarg);             break;
                case 'y':
                        u = stod(optarg); 		break;
	        case 'z':
                        g = stod(optarg); 		break;
	        case 'j':
                        lambda = stod(optarg);          break;
                case 'b':
                        kappa = stod(optarg);           break;
                case 'k':
                        lambda0 = stod(optarg);         break;
                case 'l':
                        lambdax = stod(optarg);         break;
                case 'm':
                        lambdaz = stod(optarg);         break;
                case 'i':
                        K0 = stod(optarg);              break;
                case 'n':
                        passes = stoi(optarg);          break;
                case 's':
                        if (soptarg.compare("True") == 0)
                        {
                                printLastConf_Q = "True";
                        } else {
                                printLastConf_Q = "False";
                        }				break;
                case 'r':
                        rs = stoi(optarg);	 	break;
                case '?':
                        /* getopt_long already printed an error message. */ break;
        default:
          abort ();
        }
    }
}



//Few useful functions to start with:
// random real between 0 and 1
double randomreal() {return (double) rand() / RAND_MAX;}

// random number between -1 and 1
double randomrealpm() {return ((double) 2.0 * randomreal() - 1.0);}


// get vector index from 2d array (i,j) \elem [0,size-1]^3  ->  k \elem [0, size^3 - 1] ...
int index(int i, int j) {return i + j * size;}

// ... and back
int cx(int k) {return k % size;}
int cy(int k) {return (int) (k - cx(k)) / size;}





// Define Energy of the model
// Onsite energy
double hamiltonian_onsite(double R1, double I1, double R2, double I2, int site_index)
{
    double e = 0.0;

    e += r0 * (R1 * R1 + I1 * I1 + R2 * R2 + I2 * I2) / 2;

    e += -x * (R1 * R1 + I1 * I1 - R2 * R2 - I2 * I2) / 2;

    e += eps * (R1 * R2 + I1 * I2);

    e += (u + lambda) * (R1 * R1 + I1 * I1 + R2 * R2 + I2 * I2) * (R1 * R1 + I1 * I1 + R2 * R2 + I2 * I2) / 8;

    e += (u - lambda) * (R1 * R1 + I1 * I1 - R2 * R2 - I2 * I2) * (R1 * R1 + I1 * I1 - R2 * R2 - I2 * I2) / 8;

    e += g * (R1 * R2 + I1 * I2) * (R1 * R2 + I1 * I2) / 2;

    e += kappa * (R1 * R1 + I1 * I1 + R2 * R2 + I2 * I2) * (R1 * R1 + I1 * I1 - R2 * R2 - I2 * I2) / 4;

    e += lambda0 * f[site_index] * (R1 * R1 + I1 * I1 + R2 * R2 + I2 * I2);

    e += lambdax * f[site_index] * (R1 * R2 + I1 * I2);

    e += lambdaz * f[site_index] * (R1 * R1 + I1 * I1 - R2 * R2 - I2 * I2);
    
    return e;
}

//Nearest neighbor interaction
double hamiltonian_nn(double R1, double I1, double R2, double I2, int site_index)
{
    double e = 0.0;
        
    double nR1;
    double nI1;
    double nR2;
    double nI2;
    
    nR1 = R1_state[nnup[site_index]];
    nI1 = I1_state[nnup[site_index]];
    nR2 = R2_state[nnup[site_index]];
    nI2 = I2_state[nnup[site_index]];

    //Interaction with nearest neighbor
    e += K0 * ( (R1 - nR1) * (R1 - nR1) + (I1 - nI1) * (I1 - nI1) + (R2 - nR2) * (R2 - nR2) + (I2 - nI2) * (I2 - nI2));

    nR1 = R1_state[nndown[site_index]];
    nI1 = I1_state[nndown[site_index]];
    nR2 = R2_state[nndown[site_index]];
    nI2 = I2_state[nndown[site_index]];

    //Interaction with nearest neighbor
    e += K0 * ( (R1 - nR1) * (R1 - nR1) + (I1 - nI1) * (I1 - nI1) + (R2 - nR2) * (R2 - nR2) + (I2 - nI2) * (I2 - nI2));

    nR1 = R1_state[nnleft[site_index]];
    nI1 = I1_state[nnleft[site_index]];
    nR2 = R2_state[nnleft[site_index]];
    nI2 = I2_state[nnleft[site_index]];

    //Interaction with nearest neighbor
    e += K0 * ( (R1 - nR1) * (R1 - nR1) + (I1 - nI1) * (I1 - nI1) + (R2 - nR2) * (R2 - nR2) + (I2 - nI2) * (I2 - nI2));

    nR1 = R1_state[nnright[site_index]];
    nI1 = I1_state[nnright[site_index]];
    nR2 = R2_state[nnright[site_index]];
    nI2 = I2_state[nnright[site_index]];

    //Interaction with nearest neighbor
    e += K0 * ( (R1 - nR1) * (R1 - nR1) + (I1 - nI1) * (I1 - nI1) + (R2 - nR2) * (R2 - nR2) + (I2 - nI2) * (I2 - nI2));

    // I calculated the gradient twice each time
    e = e / 2;

    return e;
}

//Full hamiltonian
double hamiltonian(double R1, double I1, double R2, double I2, int site_index)
{
    double e = 0.0;
        
    e += hamiltonian_onsite(R1, I1, R2, I2, site_index);
    
    e += hamiltonian_nn(R1, I1, R2, I2, site_index);
    
    return e;
}





// Initialize a random spin configuration, calculate its energy & magnetization
void initialization ()
{
	// set the random number generator seed
	srand(rs);
    
    // Generate an initial spin configuration
	// Evaluate for each site the neighboring site indices
	// Define sweep order (initially ordered 0, 1, 2, ...)
	for (int i = 0; i < N; i++)
	{
            double theta = acos( randomrealpm() );
	    	    
            R1_state[i] = randomreal();
            I1_state[i] = randomrealpm();
            R2_state[i] = randomrealpm();
            I2_state[i] = randomrealpm();

            f[i] = 1.0 * (cx(i) - (1.0 * size)/2 ) / ((cx(i) - (1.0 * size)/2 ) * (cx(i) - (1.0 * size)/2 ) + (cy(i) - (1.0 * size)/2 ) * (cy(i) - (1.0 * size)/2 ) + 1.0e-9);
	    
            nnup[i]     = index(cx(i), (cy(i) + 1) % size);
	    nndown[i]   = index(cx(i), (cy(i) + size - 1) % size);
            nnleft[i]   = index((cx(i) + size - 1) % size, cy(i));
	    nnright[i]  = index((cx(i) + 1) % size, cy(i));
	    sweeporder[i] = i;
	}

	// Evaluate the Energy and Magnetization of the initial state
	for (int i = 0; i < N; i++)
	{
            double R1 = R1_state[i];
            double I1 = I1_state[i];
            double R2 = R2_state[i];
            double I2 = I2_state[i];
	    	    
            E += hamiltonian_nn(R1, I1, R2, I2, i) / 2 + hamiltonian_onsite(R1, I1, R2, I2, i);
	}
	
        // All quantities are per spin
	E /= N;
}



//All necessary functions for printing out useful output
//Generate filename for the main data file
string filename()
{
        string fn = "dislocation-solver_r0=" + to_string(r0).substr(0,6) + "_lambda0=" + to_string(lambda0).substr(0,6) + "_lambdax=" + to_string(lambdax).substr(0,6) + "_lambdaz=" + to_string(lambdaz).substr(0,6) + "_K0=" + to_string(K0).substr(0,6) + "_passes=" + to_string(passes) + "_rs=" + to_string(rs);
	return fn;
}



// Create output file and write header line
void init_printf()
{
	ofstream outfile;
	string fn = filename() + ".tsv";
	
	outfile.open(fn);
        outfile << "pass \t Temperature \t x \t eps \t u \t g \t lambda \t kappa \t lambda0 \t lambdax \t lambdaz \t K0 \t <E> \t <E^2>\n";
	outfile.close();
}

// Export data to file
void printf(double e, double e2)
{
	ofstream outfile;
	string fn = filename() + ".tsv";
	
	outfile.open(fn,ios_base::app);
        outfile << to_string(pass) + "\t" + to_string(r0) + "\t" + to_string(x) + "\t" + to_string(eps) + "\t" + to_string(u) + "\t" + to_string(g) + "\t" + to_string(lambda) + "\t" + to_string(kappa) + "\t" + to_string(lambda0) + "\t" + to_string(lambdax) + "\t" + to_string(lambdaz) +  "\t" + to_string(K0) +"\t" + to_string(e) + "\t" + to_string(e2) + "\n";
	outfile.close();
}








// Create output file and write header line
void init_printlog()
{
	ofstream outfile;
	string fn = filename() + ".log";
	
	outfile.open(fn);
	outfile << "Starting the simulation with the parameters \n";
        outfile << 	"r0 = "  + to_string(r0) + "\t" +
                        "x = " + to_string(x) + "\t" +
                        "eps = " + to_string(eps) + "\t" +
                        "u = " + to_string(u) + "\t" +
                        "g = " + to_string(g) + "\t" +
                        "lambda = " + to_string(lambda) + "\t" +
                        "kappa = " + to_string(kappa) + "\t" +
                        "lambda0 = " + to_string(lambda0) + "\t" +
                        "lambdax = " + to_string(lambdax) + "\t" +
                        "lambdaz = " + to_string(lambdaz) + "\t" +
                        "K0 = " + to_string(K0) + "\t" +
                        "passes = " + to_string(passes) + "\n\n";
	outfile.close();

}

// Export data to file
void printlog(string logstr)
{
	ofstream outfile;
	string fn = filename() + ".log";
	
	outfile.open(fn,ios_base::app);	
	outfile << logstr;
	outfile.close();
}








//Generate filename for the last configuration data file
string filenameLC()
{
        string fn = filename() + "_finalstate" + ".tsv";
	return fn;
}

// Create output file for last configuration
void init_printLC()
{
	ofstream outfile;
	string fn = filenameLC();
	
	outfile.open(fn);
        outfile << "State of a " << size << "x" << size << " system  (R10, I10, R20, I20), (R11, I11, R21, I21), ... \n";
        //In Mathematica one can use Partition[list,size] to make this list square again
	outfile.close();
}

void printLC()
{
	ofstream outfile;
	string fn = filenameLC();

	outfile.open(fn,ios_base::app);
	for (int i = 0; i < N; i++)
	{
                outfile << to_string(R1_state[i]) << "\t" << to_string(I1_state[i]) << "\t" << to_string(R2_state[i]) << "\t" << to_string(I2_state[i]) << "\n";
	}
	outfile.close();
}





// All necessary functions to generate one pass:
// Generate a random sweep order to visit the sites
void generate_new_sweep_order()
{
	// Shuffle the sweep order
	for (int i = 0; i < N; i++)
	{
	    int mem = sweeporder[i];
	    int newi = random() % N;
	    
	    sweeporder[i] = sweeporder[newi];
	    sweeporder[newi] = mem;
	}
}


// define one spin flip (Metropolis-Hastings step)
// This step also accounts for possible phase space reduction
void oneflip(int i)
{
        double R1 = R1_state[i];
        double I1 = I1_state[i];
        double R2 = R2_state[i];
        double I2 = I2_state[i];

        double nR1 = R1 + opening * randomrealpm();
        double nI1 = I1 + opening * randomrealpm();
        double nR2 = R2 + opening * randomrealpm();
        double nI2 = I2 + opening * randomrealpm();

        double delta_E = hamiltonian(nR1, nI1, nR2, nI2, i) - hamiltonian(R1, I1, R2, I2, i);
	

        //watch out for the choice of MetropolisT
        double MetropolisT = 0.001 * x;
        if ((delta_E < 0) || (randomreal() < exp(- delta_E / MetropolisT))) {
        
        R1_state[i] = nR1;
        I1_state[i] = nI1;
        R2_state[i] = nR2;
        I2_state[i] = nI2;

        E += delta_E / N;
        new_acceptance_rate += 1;
    } 
}



// Defines a full pass of spin flips
// This function evaluates the new acceptance rate and adjusts  the opening angle if necessary
void onepass()
{
	generate_new_sweep_order();
	
	new_acceptance_rate = 0;
	for (int i = 0; i < N; i++)
	{
        oneflip(sweeporder[i]);
    }

    acceptance_rate = new_acceptance_rate / N;
    
    if (acceptance_rate < thr_acceptance_rate)
    {
        opening = max( opening / 1.01, 0.0);
    } else
    {
        opening = min( opening * 1.01, pi);
    }

    for (int i = 0; i < size; i++)
    {
        I1_state[index(i,0)] = 0;
        I1_state[index(i,size-1)] = 0;
        I1_state[index(0,i)] = 0;
        I1_state[index(size-1,i)] = 0;
    }
}

void compute()
{
	// roughly half of the passes should be added as prepasses
	int prepasses = (int) passes / 2;
	if (prepasses > 500000)
		prepasses = 500000;
	
        // after performing all subpasses, an output data is generated
	int subpasses = 100000;
	if (passes <= subpasses){subpasses = (int) passes/10;}
	
        // conversely one gets one output data per superpass
	int superpasses = (int) passes / subpasses;
	
	
        // compute prepasses
	pass = 0;
	for (int i = 0; i < prepasses; i++)
	{
		onepass();
		pass += 1;
	}
	

        //initialize measurable computation
	init_printf();
	init_printlog();
	
	for (int i = 0; i < superpasses; i++)
	{
		double e = 0.0;
		double e2 = 0.0;

		av_acceptance_rate = 0.0;
		for (int j = 0; j < subpasses; j++)
		{
			onepass();
			e += E;
			e2 += E * E;

			pass += 1;
			av_acceptance_rate += acceptance_rate;	
		}
		e /= subpasses;
		e2 /= subpasses;

		av_acceptance_rate /= subpasses;

                printf(e, e2);
                printlog("pass = " + to_string(pass) + "\t" + "av_acceptance_rate = " + to_string(av_acceptance_rate) + "\t" + "opening = " + to_string(opening) + "\n");
	}
	
	if (printLastConf_Q.compare("True") == 0)
	{
		init_printLC();
		printLC();
	}
}

int main(int argc, char** argv)
//-------
{
    parse_input(argc, argv);
    
    cout << "Starting the simulation with the parameters" <<  endl;
    cout << "r0 = "  + to_string(r0) + "\t" +
            "x = " + to_string(x) + "\t" +
            "eps = " + to_string(eps) + "\t" +
            "u = " + to_string(u) + "\t" +
            "g = " + to_string(g) + "\t" +
            "lambda = " + to_string(lambda) + "\t" +
            "kappa = " + to_string(kappa) + "\t" +
            "lambda0 = " + to_string(lambda0) + "\t" +
            "lambdax = " + to_string(lambdax) + "\t" +
            "lambdaz = " + to_string(lambdaz) + "\t" +
            "K0 = " + to_string(K0) + "\t" +
            "passes = " + to_string(passes) + "\n\n" << endl;
    time_t t0 = time(nullptr);
	
    initialization();
	
    compute();
	
    time_t t = time(nullptr);
    double dt = 1.0 * (t - t0);
    cout << "Simulation time = " << dt << endl;
    printlog("Overall simulation time = " + to_string(dt) + "\n");

    return 0;
}
