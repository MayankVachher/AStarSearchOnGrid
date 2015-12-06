#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <limits>
#include <queue>
#include <iostream>

/* To Calculate Memory Usage
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( ) {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}





/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( ) {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}

// Main Algorithm starts here
// Author: Mayank Vachher

using namespace std;
double inf = std::numeric_limits<double>::infinity();
double time_taken = 0;
int statesExpanded = 0;

class node { // A node class
	node *par; // Parent node

	double genHeuristic(int *loc, int g_i, int g_j) {
		double dx = abs(loc[0] - g_i);
    	double dy = abs(loc[1] - g_j);
		return dx + dy - 0.586 * min(dx,dy); // Octile Distance
	}

	public:
	int *loc; // self position
	int val; // value of node
	double f, g, h; // function values

	node* getParent() {
		return par;
	}

	node(int i, int j, int got_val) {
		val = got_val;
		loc = new int[2];
		loc[0] = i;
		loc[1] = j;
		g = inf;
		h = 0;
		f = g + h;
		par = NULL;
	}

	node(int i, int j, int got_val, int g_i, int g_j) {
		val = got_val;
		loc = new int[2];
		loc[0] = i;
		loc[1] = j;
		par = NULL;
		g = inf;
		h = genHeuristic(loc, g_i, g_j);
		f = g + h;
	}

	bool update_g(double val) { // update the g value
		if(val < g) { // only update it if you can make it better
			g = val;
			f = g + h;
			return true; // tell the function call to update the parent also
		}
		return false;
	}

	void updateParent(node* parent) {
		par = parent;
	}
};

struct CompareNode : public std::binary_function<node*, node*, bool> { // for priorityQueue
	bool operator()(const node* lhs, const node* rhs) const {
		if(lhs->f!=rhs->f)
			return lhs->f > rhs->f;
		return lhs->h < rhs->h;
	}
};

template<
    class T,
    class Container = std::vector<T>,
    class Compare = std::less<typename Container::value_type>
> class MyQueue : public std::priority_queue<T, Container, Compare>
{
public:
    typedef typename
        std::priority_queue<
        T,
        Container,
        Compare>::container_type::const_iterator const_iterator;

    bool find(const T&val) const
    {
        auto first = this->c.cbegin();
        auto last = this->c.cend();
        while (first!=last) {
            if (*first==val) return true;
            ++first;
        }
        return false;
    }
};

class base {
	int rows, cols;
	node ***grid;
	public:
	base(int r, int c) {
		rows = r;
		cols = c;
		grid = new node**[r];
		for (int i = 0; i < r; ++i) {
  			grid[i] = new node*[c];
		}
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
    			grid[i][j] = NULL;
  			}
		}
	}
	void printMeta() {
		printf("%d, %d\n", rows, cols);
	}

	void setNode(int i, int j, int new_val) {
		grid[i][j] = new node(i, j, new_val);
	}

	void setNode(int i, int j, int new_val, int g_i, int g_j) {
		grid[i][j] = new node(i, j, new_val, g_i, g_j);
	}

	node* getNode(int i, int j) {
		return grid[i][j];
	}

	vector<node*> getUnitNeighbors(node* curr) {
		vector<node*> neighbors;
		int r = curr->loc[0],c = curr->loc[1];
		if (r == 0) {
			neighbors.push_back(getNode(r+1,c));
			if(c == 0) {
				neighbors.push_back(getNode(r,c+1));
			}
			else if(c == cols-1) {
				neighbors.push_back(getNode(r,c-1));
			}
			else {
				neighbors.push_back(getNode(r,c+1));
				neighbors.push_back(getNode(r,c-1));
			}
		}
		else if (r == rows-1) {
			neighbors.push_back(getNode(r-1,c));
			if(c == 0) {
				neighbors.push_back(getNode(r,c+1));
			}
			else if(c == cols-1) {
				neighbors.push_back(getNode(r,c-1));
			}
			else {
				neighbors.push_back(getNode(r,c+1));
				neighbors.push_back(getNode(r,c-1));
			}
		}
		else {
			neighbors.push_back(getNode(r-1,c));
			neighbors.push_back(getNode(r+1,c));
			if(c == 0) {
				neighbors.push_back(getNode(r,c+1));
			}
			else if(c == cols-1) {
				neighbors.push_back(getNode(r,c-1));
			}
			else {
				neighbors.push_back(getNode(r,c+1));
				neighbors.push_back(getNode(r,c-1));
			}			
		}
		return neighbors;
	}

	vector<node*> getDiagNeighbors(node* curr) {
		vector<node*> neighbors;
		int r = curr->loc[0],c = curr->loc[1];
		if (r == 0) {
			if(c == 0) {
				neighbors.push_back(getNode(r+1,c+1));
			}
			else if(c == cols-1) {
				neighbors.push_back(getNode(r+1,c-1));
			}
			else {
				neighbors.push_back(getNode(r+1,c+1));
				neighbors.push_back(getNode(r+1,c-1));
			}
		}
		else if (r == rows-1) {
			if(c == 0) {
				neighbors.push_back(getNode(r-1,c+1));
			}
			else if(c == cols-1) {
				neighbors.push_back(getNode(r-1,c-1));
			}
			else {
				neighbors.push_back(getNode(r-1,c+1));
				neighbors.push_back(getNode(r-1,c-1));
			}
		}
		else {
			if(c == 0) {
				neighbors.push_back(getNode(r+1,c+1));
				neighbors.push_back(getNode(r-1,c+1));
			}
			else if(c == cols-1) {
				neighbors.push_back(getNode(r+1,c-1));
				neighbors.push_back(getNode(r-1,c-1));
			}
			else {
				neighbors.push_back(getNode(r+1,c+1));
				neighbors.push_back(getNode(r-1,c+1));
				neighbors.push_back(getNode(r+1,c-1));
				neighbors.push_back(getNode(r-1,c-1));
			}			
		}
		return neighbors;
	}

	double solve(node *start, node *goal) { // A* Algorithm
		if(start->val == 1 || goal->val == 1)
			return 10000000;
		MyQueue<node*, vector<node*>, CompareNode > closed_set, open_set;
		start->update_g(0);
		open_set.push(start);
		node *current_node;
		while(!open_set.empty()) {
			current_node = open_set.top(); open_set.pop();
			if(current_node == goal) {
				statesExpanded = closed_set.size() + 1;
				return current_node->f;
			}
			closed_set.push(current_node);
			vector<node*> neighbors = getUnitNeighbors(current_node);
			// printf("Curr F: %.3f\n", current_node->f);
			for(node* neighbor : neighbors) {
				// printf("Unit: %d %d\n", neighbor->loc[0], neighbor->loc[1]);
				if(closed_set.find(neighbor) || neighbor->val == 1)
					continue;
				if(neighbor->update_g(current_node->g + 1)) {
					neighbor->updateParent(current_node);
				}
				if(!open_set.find(neighbor)) {
					open_set.push(neighbor);
				}
			}

			neighbors = getDiagNeighbors(current_node);
			for(node* neighbor : neighbors) {
				// printf("Diag: %d %d\n", neighbor->loc[0], neighbor->loc[1]);
				if(closed_set.find(neighbor) || neighbor->val == 1)
					continue;
				if(neighbor->update_g(current_node->g + 1.414)) {
					neighbor->updateParent(current_node);
				}
				if(!open_set.find(neighbor)) {
					open_set.push(neighbor);
				}
			}

			// printf("\n\n");
		}
		statesExpanded = closed_set.size();
		return 10000000;

	}

	void writeSolution(node *start, node *goal) {
		FILE *toWrite = fopen("Solution.txt", "w");
		fclose(toWrite);
		if(goal->getParent() != NULL) {
			recursiveSolve(start, goal);
		}
	}
	void recursiveSolve(node* s, node* c) {
		if(c == s) {
			FILE *toWrite = fopen("Solution.txt", "a");
			fprintf(toWrite, "%d,%d\n", c->loc[1],c->loc[0]);
			fclose(toWrite);
		}
		else {
			recursiveSolve(s, c->getParent());
			FILE *toWrite = fopen("Solution.txt", "a");
			fprintf(toWrite, "%d,%d\n", c->loc[1],c->loc[0]);
			fclose(toWrite);
		}
	}
};

int main(int argc, char *argv[]) {
	if(argc < 2 || argc > 2) {
		printf("Invalid call to solver...\nUsage: ./<solver_exec> <file_to_be_solved>\n\n");
		exit(0);
	}
	int r,c, i, j, temp;
	int start_i, start_j, goal_i, goal_j;
	node *start, *goal;
	char buff[30];

	clock_t tStart = clock(); // start clock

	FILE *readFile = fopen(argv[1],"r"); // Open the file to Read
	fscanf(readFile,"%s %d %d", buff, &c, &r); // Read Grid Specs
	i = r, j = c;
	base *ENV = new base(r,c); // Create an Empty Grid
	fscanf(readFile,"%s %d %d", buff, &start_j, &start_i); // Read Start position
	fscanf(readFile,"%s %d %d", buff, &goal_j, &goal_i); // Read Goal position
	fscanf(readFile, "%s", buff); // Read Environment Line
	r = 0;
	c = 0;
	while(r<i) {
		c = 0;
		while(c<j) {
			fscanf(readFile, "%d", &temp); // Read the grid from file
			ENV->setNode(r, c, temp, goal_i, goal_j); // Populate the Grid on RAM
			c++;
		}
		r++;
	}
	double solvedCost = 10000000;
	fclose(readFile); // Done reading the file, Close it.
	if(start_i < 0 || start_j < 0 || goal_i < 0 || goal_j < 0 || start_i >= r || goal_i >= r || goal_j >= c || start_j >= c) {
		FILE *toWrite = fopen("Solution.txt", "w");
		fclose(toWrite);
	}
	else {
		start = ENV->getNode(start_i, start_j); // Point to the start cell
		goal = ENV->getNode(goal_i, goal_j); // Point to the goal cell

		solvedCost = ENV->solve(start, goal); // Solve the Grid

		time_taken = (double)(clock() - tStart)/CLOCKS_PER_SEC; // stop clock
		

		ENV->writeSolution(start, goal); // Write down the solution
	}

	//print cost, number of states expanded, Time for Running the Algo, Peak memory usage of the Whole Program
	printf("Path Cost: %.3f, States Expanded: %d, Run Time: %.2f, Memory: ", solvedCost, statesExpanded, time_taken);
	printf("%.2f\n", getPeakRSS()/1000000.0);
	return 0;
}
