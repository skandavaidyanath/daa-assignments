#include <vector>
#include <bits/stdc++.h>
#include "point.h"

using namespace std;


vector<Point> read_points(string filepath)
/**
 * reads points into a vector
 */
{
    vector<Point> vec;
    ifstream infile(filepath);
    double a, b;
    while (infile >> a >> b)
    {
        vec.push_back(Point(a, b));
    }
    infile.close();
    // sort(vec.begin(), vec.end());
    return vec;
}


void write_points(vector<int> vec, string filepath)
/**
 * writes points into file
 */
{
    ofstream outfile(filepath);
    if (outfile.is_open())
    {
        for (auto i: vec)
        {
            outfile << i << "\n";
        }
        outfile.close();
    }
    else
        cerr << "Unable to open file\n";
}


double single_error(Point p, Point line1, Point line2)
/**
 * error for a point given the regression line's endpoints
 */
{
    double m = (line1.y() - line2.y()) / (line1.x() - line2.x());
    double c = line1.y() - m * line1.x();
    double pred = (m * p.x() + c);
    return (pred-p.y())*(pred-p.y());
}


double total_error(vector<Point> points, int i, int j)
/**
 * find total error for all points between i and j endpoints of regression line
 */
{
    double err = 0;
    for (int k = i; k < j+1; k++)
    {
        err += single_error(points[k], points[i], points[j]);
    }
    return err;
}


double fill_OPT(vector<Point> points, int j, double opt[], double c)
/**
 * fill DP table
 */
{
    if (opt[j] >= 0)
        return opt[j];
    double opt_j = INF;
    for(int i = 0; i < j; i++)
    {
        double temp = fill_OPT(points, i, opt, c) + c + total_error(points, i, j);
        if (temp < opt_j)
            opt_j = temp;
    }
    opt[j] = opt_j;
    return opt[j];
}


vector<int> get_lines(vector<Point> points, int j, double opt[], double c)
/**
 * given DP opt table, find the optimum endpoints
 */
{
    if (j == 0)
    {
        vector<int> vec;
        vec.push_back(j);
        return vec;
    }
    int minInd = 0;
    for(int i = 1; i < j; i++)
    {
        double temp = opt[i] + c + total_error(points, i, j);
        if (temp < opt[minInd] + c + total_error(points, minInd, j))
            minInd = i;
    }
    vector<int> vec = get_lines(points, minInd, opt, c);
    vec.push_back(j);
    return vec;
}


vector<int> get_regression_points(vector<Point> points, double c)
/**
 * function to get optimum indices for regression upon a vector of points given cost c
 */
{
    double opt[points.size()];
    opt[0] = 0;
    for(int i = 1; i < points.size(); i++)
        opt[i] = -1;
    fill_OPT(points, points.size()-1, opt, c);
    /*
    cout << "filled opt\n";
    for (int i = 0; i < points.size(); i++)
        cout << opt[i] << ' ';
    cout << "\n";
    */
    vector<int> line = get_lines(points, points.size()-1, opt, c);
    return line;
}


int main(int argc, char* argv[])
{
    double c = stod(argv[1]);
    vector<Point> points = read_points(argv[2]);
    vector<int> line = get_regression_points(points, c);
    /*
    for (auto p: line)
        cout << p << "\n";
    */
    write_points(line, "out.txt");
}
