#include <string>
#include <sstream>

using namespace std;
#define INF (numeric_limits<double>::infinity())

class Point
{
  private:
    double _x, _y;

  public:
    Point()
    {
        this->_x = 0;
        this->_y = 0;
    }
    Point(double x1, double y1)
    {
        this->_x = x1;
        this->_y = y1;
    }
    double x() const
	/**
	 * getter method for a point's x-coordinate.
	 */
    {
        return _x;
    }
    double y() const
	/**
	 * getter method for a point's y-coordinate.
	 */
    {
        return _y;
    }
    double slope(Point* p)
	/**
	 * method to calculate the slope between the point and another.
	 */
    {
        double s;
        double y_diff = this->_y - p->y();
        double x_diff = this->_x - p->x();
        if (x_diff == 0)
            s = INF;
        else
            s = y_diff / x_diff;
        return s;
    }
    double distance_sq(Point* p)
	/**
	 * method to calculate the distance squared between the point and another.
	 */
    {
        return (this->_x-p->x())*(this->_x-p->x()) + (this->_y-p->y())*(this->_y-p->y());
    }
    string to_string()
	/**
	 * returns a string representation of a point.
	 */
    {
        ostringstream strs;
        strs << _x << " " << _y;
        return strs.str();
    }
	void update(double x, double y)
	/**
	 * method to update coordinates of a point.
	 */
	{
		this->_x=x;
		this->_y=y;
	}
    static int orientation(Point* p, Point* q, Point* r)
	/**
	 * static method to calculate the orientation of three points.
	 */
    {
        double val = (q->y() - p->y()) * (r->x() - q->x()) -
                     (q->x() - p->x()) * (r->y() - q->y());

        if (val == 0)
            return 0;
        return (val > 0) ? 1 : 2;
    }
};
