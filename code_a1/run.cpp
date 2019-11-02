#include <ctime>
#include "graph.h"

int main()
{
    vector<Point *> v = Graph::read_points("TestCases/points4.txt");
    clock_t begin, end;
    double t;
    vector<Point *> hull;

    begin = clock();
    GrahamGraph gg = GrahamGraph(v);
    hull = gg.calculate_hull();
    end = clock();
    t = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Graham's Scan: " << t << "s\n";
    gg.save_points(hull, "TestCases/outGraham.txt");

    begin = clock();
    JarvisGraph jg = JarvisGraph(v);
    hull = jg.calculate_hull();
    end = clock();
    t = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Jarvis' March: " << t << "s\n";
    jg.save_points(hull, "TestCases/outJarvis.txt");

    begin = clock();
    KPGraph kpg = KPGraph(v);
    hull = kpg.calculate_hull();
    end = clock();
    t = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Kirckpatrick Siedel: " << t << "s\n";
    kpg.save_points(hull, "TestCases/outKS.txt");

    return 0;
}
