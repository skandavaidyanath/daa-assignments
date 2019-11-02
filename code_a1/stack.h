#include <vector>
#include "point.h"


class point_stack
{
  private:
    vector<Point*> stack;
  public:
    Point* pop()
    {
        Point* p = stack.back();
        stack.pop_back();
        return p;
    }
    void push(Point* p)
    {
        stack.push_back(p);
    }
    Point* top()
    {
        return stack.back();
    }
    Point* under_top()
    {
        return stack.rbegin()[1];
    }
    bool empty()
    {
        return (stack.size() == 0);
    }
};
