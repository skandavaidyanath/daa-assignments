#include <cmath>
#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include "stack.h"

class Graph
/**
 * Parent Class for all Hull algorithms.
 */
{
  protected:
	vector<Point *> points;
	/**
	 * The vector of points representing the input to the any convex hull algorithm.
	 */

  public:
	Graph(){};
	virtual void insert_point(Point *p) {}
	/**
	 * Virtual method (to be implemented by subclasses) to insert another point in the graph.
	 */

	virtual void sort_points() {}
	/**
	 * Virtual method (to be implemented by subclasses) to sort the vector of points.
	 */

	virtual vector<Point *> calculate_hull() {}
	/**
	 * Virtual method (to be implemented by subclasses) to return a vector of points,
	 * which represent the convex hull of the graph.
	 */

	static vector<Point *> read_points(string filepath)
	/**
	 * Static method for easy reading of points from a file.
	 * The file must have one point on each line.
	 * A point must be represented by an x coordinate followed by a space
	 * followed by a y coordinate.
	 */
	{
		vector<Point *> vec;
		ifstream infile(filepath);
		double a, b;
		while (infile >> a >> b)
		{
			vec.push_back(new Point(a, b));
		}
		infile.close();
		return vec;
	}

	static void save_points(vector<Point *> vec, string filepath)
	/**
	 * Static method for easy saving of points to a file.
	 * Will write points to a file line by line.
	 * A point is represented by an x coordinate followed by a space
	 * followed by a y coordinate.
	 */
	{
		ofstream outfile(filepath);
		if (outfile.is_open())
		{
			string str;

			for (auto p : vec)
			{
				outfile << p->to_string() << "\n";
			}
			outfile.close();
		}
		else
			cerr << "Unable to open file";
	}

	void print_points()
	/**
	 * Method to print the graph's points to stdout.
	 */
	{
		for (auto p : points)
		{
			cout << p->to_string() << "\n";
		}
	}
};

class GrahamGraph : public Graph
/**
 * Class representing a graph upon which the Graham's Scan algorithm
 * will be run.
 */
{
  public:
	GrahamGraph();

	GrahamGraph(vector<Point *> points1) : Graph()
	{
		points = points1;
		sort_points();
	}

	void sort_points()
	/**
	 * Sorts the points with respect to the point at the bottom left corner.
	 */
	{
		long long minInd = 0;
		for (long long i = 1; i < points.size(); i++)
		{
			Point *min = points[minInd];
			double y = points[i]->y();
			if ((y < min->y()) || (min->y() == y && points[i]->x() < min->x()))
				minInd = i;
		}
		iter_swap(points.begin(), points.begin() + minInd);

		sort(points.begin() + 1, points.end(), [this](const auto &lhs, const auto &rhs) {
			int o = Point::orientation(points[0], lhs, rhs);
			if (o == 0)
				return (points[0]->distance_sq(rhs) >= points[0]->distance_sq(lhs));
			return (o == 2);
		});
	}

	void insert_point(Point *p)
	/**
	 * Insert a point, preservint the order of points w.r.t. the bottom left point.
	 */
	{
		if (points.size() == 0)
			points.push_back(p);
		return;
		Point *min_p = points[0];
		int left = 0, right = points.size() - 1, ind = 0;
		while (left <= right)
		{
			ind = left + (right - left) / 2;
			if (points[ind]->slope(min_p) == p->slope(min_p))
				break;
			if (p->slope(min_p) > points[ind]->slope(min_p))
				left = ind + 1;
			else
				right = ind - 1;
		}
		while (ind > 0 && p->slope(min_p) <= points[ind]->x())
			ind--;
		while (ind < points.size() && p->slope(min_p) > points[ind]->slope(min_p))
			ind++;
		points.insert(points.begin() + ind, p);
	}

	vector<Point *> calculate_hull()
	{
		vector<Point *> new_points;
		new_points.push_back(points[0]);
		for (long long i = 1; i < points.size(); i++)
		{
			while (i < points.size() - 1 && Point::orientation(points[0], points[i], points[i + 1]) == 0)
				i++;

			new_points.push_back(points[i]);
		}

		point_stack S;
		S.push(new_points[0]);
		S.push(new_points[1]);
		S.push(new_points[2]);

		for (long long i = 3; i < new_points.size(); i++)
		{
			while (Point::orientation(S.under_top(), S.top(), new_points[i]) != 2)
				S.pop();
			S.push(new_points[i]);
		}

		vector<Point *> hull;
		while (!S.empty())
		{
			hull.push_back(S.pop());
		}
		return hull;
	}
};

class JarvisGraph : public Graph
/**
 * Class representing a graph upon which the Jarvis March algorithm
 * will be run.
 */
{
  protected:
	Point *leftmost;
	void update_leftmost()
	/**
	 * update the leftmost point from the entire vector of points.
	 */
	{
		for (auto p : points)
		{
			if (p->x() < leftmost->x())
				leftmost = p;
		}
	}

	void update_leftmost(Point *p)
	/**
	 * depending upon the point given, update the leftmost point.
	 */
	{
		if (p->x() < leftmost->x())
			leftmost = p;
	}

  public:
	JarvisGraph();

	JarvisGraph(vector<Point *> points1) : Graph()
	{
		points = points1;
		leftmost = new Point(INF, INF);
		update_leftmost();
	}

	void insert_point(Point *p)
	/**
	 * Insert a point, updating the leftmost point if needed.
	 */
	{
		points.push_back(p);
		update_leftmost(p);
	}

	vector<Point *> calculate_hull()
	{
		vector<Point *> hull;
		if (points.size() < 3)
			return hull;

		long long l = 0;
		for (long long i = 1; i < points.size(); i++)
			if (points[i]->x() < points[l]->x())
				l = i;

		long long p = l, q;
		do
		{
			hull.push_back(points[p]);
			q = (p + 1) % points.size();
			for (long long i = 0; i < points.size(); i++)
			{
				if (Point::orientation(points[p], points[i], points[q]) == 2)
					q = i;
			}
			p = q;

		} while (p != l);
		return hull;
	}
};

class KPGraph : public Graph
/**
 * Class representing a graph upon which the Kirckpatrick Siedel algorithm
 * will be run.
 */
{
  private:
	static bool cmpfunc(Point *a, Point *b)
	/**
	 * comparer function to sort vector of points.
	 */
	{
		return (a->x() < b->x());
	}

	double findMedian(double arr[], long long n)
	/**
	 * find the median, given n
	 */
	{
		std::sort(arr, arr + n);
		return arr[n / 2];
	}

	void swap(double *a, double *b)
	/**
	 * swap two variables' values.
	 */
	{
		double temp = *a;
		*a = *b;
		*b = temp;
	}

	int partition(double arr[], long long l, long long r, long long x)
	/**
	 * partition the array.
	 */
	{
		long long i;
		for (i = l; i < r; i++)
			if (arr[i] == x)
				break;
		swap(&arr[i], &arr[r]);
		i = l;
		for (long long j = l; j <= r - 1; j++)
		{
			if (arr[j] <= x)
			{
				swap(&arr[i], &arr[j]);
				i++;
			}
		}
		swap(&arr[i], &arr[r]);
		return i;
	}

	int kthSmallest(double arr[], long long l, long long r, long long k)
	/**
	 * Find the k'th smallest value in the array. Used to compute Median.
	 */
	{
		if (k > 0 && k <= r - l + 1)
		{
			long long n = r - l + 1;
			long long i;
			double median[(n + 4) / 5];
			for (i = 0; i < n / 5; i++)
				median[i] = findMedian(arr + l + i * 5, 5);
			if (i * 5 < n)
			{
				median[i] = findMedian(arr + l + i * 5, n % 5);
				i++;
			}
			double medOfMed = (i == 1) ? median[i - 1] : kthSmallest(median, 0, i - 1, i / 2);
			long long pos = partition(arr, l, r, medOfMed);
			if (pos - l == k - 1)
				return arr[pos];
			if (pos - l > k - 1)
				return kthSmallest(arr, l, pos - 1, k);
			return kthSmallest(arr, pos + 1, r, k - pos + l - 1);
		}
		return -1.0;
	}

	bool double_equals(double a, double b, double epsilon = 0.0000001)
	/** 
	 * Compare two doubles of decimal precision 7
	 */
	{
		return std::abs(a - b) < epsilon;
	}

	double point_area(Point *p1, Point *p2, Point *p3)
	/** 
	 * Computes the area formed by the points
	 */
	{
		double t1 = p1->x() * (p2->y() - p3->y());
		double t2 = p1->y() * (p2->x() - p3->x());
		double t3 = (p2->x() * p3->y()) - (p2->y() * p3->x());
		return t1 - t2 + t3;
	}

	std::vector<Point *> UpperHull(Point *PMin, Point *PMax, std::vector<Point *> T)
	/** 
	 * Recursive function that computes the list of upper hull points
	 */
	{
		if (double_equals(PMin->x(), PMax->x()) && double_equals(PMin->y(), PMax->y())) // The base case where Pmin is returned if Pmin = Pmax
		{
			std::vector<Point *> base;
			base.push_back(PMin);
			return base;
		}
		else
		{
			std::sort(T.begin(), T.end(), cmpfunc); //sort the Points by increasing order of x coordinates
			long long len = T.size(), index;
			double a;
			// Computing the median of the x-coordinates of the points
			if (len % 2 == 0)
			{
				a = (double)(T[len / 2 - 1]->x() + T[len / 2]->x()) / 2.0;
				index = len / 2 - 1;
			}
			else
			{
				a = (T[len / 2]->x());
				index = len / 2;
			}
			std::vector<Point *> Tleft, Tright;
			std::vector<Point *> Bridge = UpperBridge(T, a); //Passing the list of points and the median of x-coordiantes of points as arguments to Upper Bridge function
			Point *Pl = Bridge[0];
			Point *Pr = Bridge[1];
			Tleft.push_back(Pl);
			Tright.push_back(Pr);
			Tleft.push_back(PMin);
			Tright.push_back(PMax);
			for (long long i = 0; i < len; ++i)
			{
				if (point_area(PMin, Pl, T[i]) > 0) // Pushing the point to Tleft if it is to the left of the line joining PMin and Pl
					Tleft.push_back(T[i]);
				if (point_area(PMax, Pr, T[i]) < 0) // Pushing the point to Tright if it is to the right of the line joining PMax and Pr
					Tright.push_back(T[i]);
			}

			std::vector<Point *> U1 = UpperHull(PMin, Pl, Tleft);  //Upper Hull on the points to on and to the left of the line joining PMin and Pl
			std::vector<Point *> U2 = UpperHull(Pr, PMax, Tright); //Upper Hull on the points to on and to the right of the line joining Pr and PMax
			//Appends the two lists U1 and U2 and returns
			for (long long i = 0; i < U2.size(); ++i)
			{
				U1.push_back(U2[i]);
			}
			return U1;
		}
	}

	std::vector<Point *> UpperBridge(std::vector<Point *> T, double a)
	/**
	 * Recursive function that computes the upper bridge of the list of points
	 */
	{
		try
		{
			std::vector<Point *> candidates;
			std::sort(T.begin(), T.end(), cmpfunc); //Sort the points in the increasing order of x-coordinates
			long long len = T.size();
			if (len == 2) // If there are only two points the line joining them is the upper bridge
			{
				std::vector<Point *> base;
				base.push_back(T[0]);
				base.push_back(T[1]);
				return base;
			}
			else
			{
				long long ind = 0;
				if (len % 2 != 0) //If there are odd number of points,push the unpaired point to candidates list
				{
					candidates.push_back(T[0]);
					ind = 1;
				}
				std::vector<pair<Point *, Point *>> Pairs;
				for (long long i = ind; i < len; i += 2)
				{
					Pairs.push_back(make_pair(T[i], T[i + 1])); //Pairing the two consecutive points
				}
				double k[Pairs.size() + 3], kk[Pairs.size() + 3];
				long long klen = 0;
				for (long long i = 0; i < Pairs.size(); i++)
				{
					if (double_equals(Pairs[i].first->x(), Pairs[i].second->x()))
					{
						//If the x-coordinates of points of a pair are equal then push the point with greater Y coordinate
						if (Pairs[i].first->y() > Pairs[i].second->y())
						{
							candidates.push_back(Pairs[i].first);
						}
						else
						{
							candidates.push_back(Pairs[i].second);
						}
					}
					else
					{
						//If slopes are not infinite compute the slope of the points of a pair and store them
						k[i] = (Pairs[i].first->y() - Pairs[i].second->y()) / (Pairs[i].first->x() - Pairs[i].second->x());
						kk[i] = k[i];
						klen++;
					}
				}
				double kmedian;
				std::sort(kk, kk + klen);
				//Finding the median of the slopes of points of pairs
				if (klen % 2 == 0)
				{
					//kmedian = kthSmallest(kk, 0, klen-1, klen/2 - 1);
					kmedian = (kk[klen / 2 - 1] + kk[klen / 2]) / 2;
				}
				else
				{
					//kmedian = kthSmallest(kk, 0, klen-1, klen/2);
					kmedian = (kk[klen / 2]);
				}
				std::vector<pair<Point *, Point *>> Small;
				std::vector<pair<Point *, Point *>> Equal;
				std::vector<pair<Point *, Point *>> Large;
				for (long long i = 0; i < klen; i++)
				{
					if (k[i] < kmedian)
						Small.push_back(Pairs[i]); //storing the pairs with slopes less than median in Small
					else if (double_equals(k[i], kmedian))
						Equal.push_back(Pairs[i]); //storing the pairs with slopes equal than median in Equal
					else
						Large.push_back(Pairs[i]); //storing the pairs with slopes greater than median in Large
				}
				double max_intercept = LLONG_MIN, intercept[T.size() + 2];
				//Finding the points with maximum intercept
				for (long long i = 0; i < T.size(); i++)
				{
					intercept[i] = T[i]->y() - kmedian * T[i]->x();
					max_intercept = max(max_intercept, intercept[i]);
				}

				Point *Pk;
				Point *Pm;
				double maxpm = LLONG_MIN, minpk = LLONG_MAX;
				for (long long i = 0; i < T.size(); ++i)
				{
					if (double_equals(intercept[i], max_intercept))
					{
						if (maxpm < T[i]->x()) // Among the points with max-intercept finding the points with maximum x-coordinate
						{
							maxpm = T[i]->x();
							Pm = T[i];
						}
						if (minpk > T[i]->x()) // Among the points with max-intercept finding the points with minimum x-coordinate
						{
							minpk = T[i]->x();
							Pk = T[i];
						}
					}
				}
				//If x(pk) < a and x(pm) > a, then return (pk, pm)
				if ((Pk->x() < a || double_equals(Pk->x(), a)) && Pm->x() > a)
				{
					std::vector<Point *> base;
					base.push_back(Pk);
					base.push_back(Pm);
					return base;
				}
				/*If x(pm)<= a then
				for all (pi, pj) belongs to  Large U Equal insert pj into candidates
				for all (pi, pj) belongs to Small insert pi and pj into candidates
				*/
				if ((Pm->x() < a || double_equals(Pm->x(), a)))
				{
					for (long long i = 0; i < Large.size(); ++i)
					{
						candidates.push_back(Large[i].second);
					}

					for (long long i = 0; i < Equal.size(); ++i)
					{
						candidates.push_back(Equal[i].second);
					}
					for (long long i = 0; i < Small.size(); i++)
					{
						candidates.push_back(Small[i].first);
						candidates.push_back(Small[i].second);
					}
				}
				/*If x(pk)> a then
				for all (pi, pj) belongs to  Small U Equal insert pi into candidates
				for all (pi, pj) belongs to Large insert pi and pj into candidates
				*/
				if (Pk->x() > a)
				{
					for (long long i = 0; i < Small.size(); ++i)
					{
						candidates.push_back(Small[i].first);
					}
					for (long long i = 0; i < Equal.size(); ++i)
					{
						candidates.push_back(Equal[i].first);
					}
					for (long long i = 0; i < Large.size(); i++)
					{
						candidates.push_back(Large[i].first);
						candidates.push_back(Large[i].second);
					}
				}

				return (UpperBridge(candidates, a)); //Returning by passing candidates list and median of x-coordinates as arguments
			}
		}
		catch (std::exception &e)
		{
			std::cout << e.what() << "\n";
		}
	}

	std::vector<Point *> solveUpper(std::vector<Point *> points)
	/**
	 * Recursive function which finds the points in the upper hull of given set of points by calling Upper_Hull
	 */
	{
		std::vector<Point *> PMin;
		std::vector<Point *> PMax;
		Point *PUMin;
		Point *PUMax;
		double minx = LLONG_MAX, maxx = LLONG_MIN, maxy1 = LLONG_MIN, maxy2 = LLONG_MIN, miny1 = LLONG_MAX, miny2 = LLONG_MAX;
		for (long long i = 0; i < points.size(); ++i)
		{
			minx = min(minx, points[i]->x());
			maxx = max(maxx, points[i]->x());
		}
		for (long long i = 0; i < points.size(); ++i)
		{
			//PUMin is the point with maximum y-coordiante among the points with minimum x-coordinate
			if (double_equals(points[i]->x(), minx))
			{
				PMin.push_back(points[i]);
				if (points[i]->y() > maxy1)
				{
					maxy1 = points[i]->y();
					PUMin = points[i];
				}
			}
			//PUMax is the point with maximum y-coordiante among the points with maximum x-coordinate
			if (double_equals(points[i]->x(), maxx))
			{
				PMax.push_back(points[i]);
				if (points[i]->y() > maxy2)
				{
					maxy2 = points[i]->y();
					PUMax = points[i];
				}
			}
		}
		std::vector<Point *> T;
		T.push_back(PUMax);
		T.push_back(PUMin);
		for (int i = 0; i < points.size(); ++i)
		{
			if (double_equals(points[i]->x(), PUMin->x())) // We only want to include points in between
			{
				continue;
			}
			if (double_equals(points[i]->x(), PUMax->x()))
			{
				continue;
			}
			T.push_back(points[i]);
		}
		std::vector<Point *> Upper = UpperHull(PUMin, PUMax, T);
		return Upper; //returns list of upper hull points
	}

  public:
	KPGraph(){};

	KPGraph(vector<Point *> points1) : Graph()
	{
		points = points1;
	}

	std::vector<Point *> calculate_hull()
	{
		if (points.size() <= 2)
		{
			std::cout << "Sorry Not sufficient no of points\n";
		}
		double minx = LLONG_MAX, miny = LLONG_MAX, maxx = LLONG_MIN, maxy = LLONG_MIN; //Finding minimum and maximum x and y coordinates
		for (long long i = 0; i < points.size(); i++)
		{
			minx = min(minx, points[i]->x());
			miny = min(miny, points[i]->y());
			maxx = max(maxx, points[i]->x());
			maxy = max(maxy, points[i]->y());
		}
		ofstream myfile;
		myfile.open("bounds.txt");
		myfile << minx - 10 << " " << miny - 10 << " " << maxx + 10 << " " << maxy + 10; // Writing the minimum and maximum x and y coordinates to file minmax.txt to plot graph
		std::vector<Point *> Hull;
		std::vector<Point *> Upper = solveUpper(points); //Computing the upper hull of the points

		//Negating the points to find the lower hull
		for (long long i = 0; i < points.size(); ++i)
		{
			points[i]->update(points[i]->x(), -1 * points[i]->y());
		}
		std::vector<Point *> Lower = solveUpper(points); //Upper hull of negated points which is same as the lower hull of points
		for (long long i = 0; i < points.size(); ++i)
		{
			points[i]->update(points[i]->x(), -1 * points[i]->y()); //negating the points back to the original
		}
		std::reverse(Lower.begin(), Lower.end());
		Hull.reserve(Upper.size() + Lower.size());
		Hull.push_back(Upper[0]);
		long long k = 0;

		//Merging the upper hull and lower hull to get the final convex hull
		for (long long i = 1; i < Upper.size(); ++i)
		{
			if (double_equals(Upper[i]->x(), Hull[k]->x()) && double_equals(Upper[i]->y(), Hull[k]->y()))
				continue;
			else
			{
				Hull.push_back(Upper[i]);
				k++;
			}
		}
		for (long long i = 0; i < Lower.size(); ++i)
		{
			if (double_equals(Lower[i]->x(), Hull[k]->x()) && double_equals(Lower[i]->y(), Hull[k]->y()))
				continue;
			else
			{
				Hull.push_back(Lower[i]); //Appending the lower hull to the end of the upper hull
				k++;
			}
		}
		if (double_equals(Hull[0]->x(), Hull[k]->x()) && double_equals(Hull[0]->y(), Hull[k]->y()))
		{
			Hull.erase(Hull.begin() + k);
		}
		return Hull;
	}
};
