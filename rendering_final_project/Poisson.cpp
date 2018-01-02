
#include <iostream>
#include <vector>
#include <random> // for random
#include <time.h> // for random
#include <fstream> // for writing results into txt file


///////////////
// CONSTANTS //
///////////////

const int NumPoints = 20000;	// minimal number of points to generate
const bool Circle = true;	// 'true' to fill a circle, 'false' to fill a rectangle
const int k = 30;		// according to bridson-siggraph07-poissondisk.pdf
const float MinDistance = sqrt(float(NumPoints)) / float(NumPoints);


//////////////////////////
// GLOBALS & STRUCTURES //
//////////////////////////

std::random_device rd;
// A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
std::mt19937 gen(rd());
// produces random floating-point value uniformly distributed on the interval [a, b).
std::uniform_real_distribution<> dis(0.0, 1.0);

float* g_DensityMap = NULL;

struct MyPoint {
	// constructor
	MyPoint() {
		x = 0; y = 0; IsValid = false;
	}
	// second constructor
	MyPoint(float a, float b) {
		x = a; y = b; IsValid = true;
	}
	// data
	float x; float y; bool IsValid;
	// methods
	bool IsInRange() const {
		bool ret = (Circle)? (((x - 0.5f) * (x - 0.5f) + (y - 0.5f) * (y - 0.5f)) <= 0.25f) : (x >= 0 && y >= 0 && x <= 1 && y <= 1);
		return ret;
	}
	// actually the following two methods are never used because we have IsInRange() already.
	bool IsInRectangle() const {
		return (x >= 0 && y >= 0 && x <= 1 && y <= 1);
	}
	bool IsInCircle() const {
		return (((x - 0.5f) * (x - 0.5f) + (y - 0.5f) * (y - 0.5f)) <= 0.25f);
	}
};

struct MyGridPoint {
	// constructor
	MyGridPoint(int a, int b) {
		x = a; y = b;
	}
	// data
	int x; int y;
};



///////////////////////
// UTILITY FUNCTIONS //
///////////////////////

/**
 * to generate a random floating-point number
 * @return a float, randomly generated.
 */
float GenerateRandomFloat() {
	return static_cast<float>(dis(gen));
}

/**
 * to get the distance between two MyPoint points
 * @param p1 the first point
 * @param p2 the second point
 * @return a floating-point number indicating the distance.
 */
float GetDistance(const MyPoint& p1, const MyPoint& p2) {
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

/**
 * to convert a MyPoint point to MyGridPoint point
 * @param p the MyPoin to be converted
 * @param CellSize the size of a cell in the space
 * @return the MyGridPoint point resulted.
 */
MyGridPoint ImageToGrid(const MyPoint& p, float CellSize) {
	return MyGridPoint((int)(p.x / CellSize), (int)(p.y / CellSize));
}

///////////////////////////
// ANOTHER STRUCTURE DEF //
///////////////////////////
// this structure is defined here because of need of ImageToGrid() utility function.
struct MyGrid {
	// constructor
	MyGrid(int _w, int _h, float _cellsize) {
		W = _w;
		H = _h;
		CellSize = _cellsize;
		Grid.resize(H);
		for(auto i=Grid.begin(); i != Grid.end(); i++) {
			i->resize(W);
		}
	}
	// methods
	void Insert(const MyPoint& p) {
		// convert to MyGridPoint from MyPoint first
		MyGridPoint gp = ImageToGrid(p, CellSize);
		Grid[gp.x][gp.y] = p;
	}
	bool IsInNeighborhood(MyPoint p, float min_dist, float CellSize)
	{
		MyGridPoint G = ImageToGrid( p, CellSize );

		// number of adjucent cells (neighboring cells) to check
		int D = 5;

		// scan the neighborhood of the point in the grid
		for(int i = (G.x - D); i < (G.x + D); i++) {
			for(int j = (G.y - D); j < (G.y + D); j++) {
				if(i >= 0 && i < W && j >= 0 && j < H) {
					// this if-statement makes sure inside the grid
					// get the point
					MyPoint p2 = Grid[i][j];
					// check if the point exists and distance is less than min_dist
					if(p2.IsValid && GetDistance(p, p2) < min_dist)
						return true;
				} else {
					// not inside the grid
					continue;
				}
			}
		}

		return false;
	}
	// data
	int W;
	int H;
	float CellSize;
	std::vector<std::vector<MyPoint>> Grid;
};


///////////////////////
// UTILITY FUNCTIONS //
// (continued)       //
///////////////////////

/**
 * randomly choose an MyPoint in the vector and pop it out
 * @param Points the vector of MyPoint
 * @return an MyPoint randomly picked and popped from the vector
 */
MyPoint PopRandom(std::vector<MyPoint>& Points) {
	std::uniform_int_distribution<> dis(0, Points.size() - 1);
	int index = dis(gen);
	MyPoint ret = Points[ index ];
	Points.erase(Points.begin() + index);
	return ret;
}

/**
 * to generate a random point that is at least MinDist from P
 * @param P the MyPoint to generate random point around
 * @param MinDist the minimum distance between the generated point and P
 * @return an MyPoint generated
 */
MyPoint GenerateRandomPointAround(const MyPoint& P, float MinDist) {
	// my own function
	float r1 = GenerateRandomFloat();
	float r2 = GenerateRandomFloat();

	// determine the radius
	// (should be between 1 * MinDist and 2 * MinDist)
	float radius = MinDist * (r1 + 1.0f);

	// determine the angle
	float theta = 2 * 3.141592653589f * r2;

	// x y coordinates for the generated point around the point (x, y)
	float x = P.x + radius * cos(theta);
	float y = P.y + radius * sin(theta);

	return MyPoint(x, y);
}

/**
 * the main function to generate samples
 * @param minDistance the minimum distance, actually derived from NumPoints (a predefined constant).
 * @param NewPointsCount how many new points (maximum) to generate (actually k = 30 in paper).
 * @param NumPoints a predefined constant, maximum count of sample points.
 * @return a vector of MyPoint with sample points inside.
 */
std::vector<MyPoint> GeneratePoissonPoints(float minDistance, int NewPointsCount, int NumPoints) {
	// the vector for processing
	std::vector<MyPoint> SamplePoints; // storing the result
	std::vector<MyPoint> ActiveList; // serve as a queue (described in the paper).

	// create the grid
	float GridCellSize = minDistance / sqrt(2.0f);
	int GridW = (int)ceil(1.0f / GridCellSize);
	int GridH = (int)ceil(1.0f / GridCellSize);

	MyGrid my_Grid(GridW, GridH, GridCellSize);

	// find the first point (generate until accepted)
	MyPoint first_point = MyPoint(GenerateRandomFloat(), GenerateRandomFloat());
	while(!(first_point.IsInRange())) {
		first_point = MyPoint(GenerateRandomFloat(), GenerateRandomFloat());
	}

	// push the first point into the vectors
	ActiveList.push_back(first_point);
	SamplePoints.push_back(first_point);
	my_Grid.Insert(first_point);

	// while queue is not empty, generate new points for each point in the queue
	while (!ActiveList.empty() && SamplePoints.size() < NumPoints) {
		// pop one point from the ActiveList
		MyPoint Point = PopRandom(ActiveList);
		// find new points
		for(int i=0; i<NewPointsCount; i++) {
			// generate one randomly
			MyPoint a_new_point = GenerateRandomPointAround(Point, minDistance);
			// check the criteria
			if(a_new_point.IsInRange() && !my_Grid.IsInNeighborhood(a_new_point, minDistance, GridCellSize)) {
				// add it (push into vectors)
				ActiveList.push_back(a_new_point);
				SamplePoints.push_back(a_new_point);
				my_Grid.Insert(a_new_point);
			}
		}
	}
	// return the resulting sample points
	return SamplePoints;
}


///////////////////
// PROGRAM ENTRY //
///////////////////

int main(int argc, char** argv) {

	// prepare my pseudorandom number generator (PRNG)
	gen.seed( time( NULL ) );

	// generate the points
	std::vector<MyPoint> Points = GeneratePoissonPoints(MinDistance, k, NumPoints);

	// dump points to a text file
	// create output file stream (std::ofstream)
	std::ofstream File("Poisson.txt", std::ios::out);

	// output the number of points
	File << Points.size() << std::endl;

	// output the coordinates of the points
	for(const auto& p : Points)
		File << "(" << p.x << ", " << p.y << ")" << std::endl;

	return 0;
}
