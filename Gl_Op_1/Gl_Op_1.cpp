//Just like in my other projects (?) there will be a lot of useless comments. Please don't mind that (As bad quality of my english as well, huh...)

#include <iostream>
#include <functional> //for std::function
#include <cmath>  //for math functions e.g. sin() or cos()
#include <vector> 
#include <algorithm>

//About style of naming. There is no style.

using std::vector;
using std::pair;

typedef pair<double, double> dpair; 



class Minimazer { //Don't ask me why... But only because using a sledge-hammer to crack a nut sounds fun.

protected:
    std::function<double(double)> func; //Function that needs findin' minimum

    double a; //Beginning of the segment
    double b; //End of the segment
    
    double r; //Coefficient of the method

    double eps;

    dpair sol;
    bool solved=false;
    unsigned long long count = 0; // How many times func was executed

public:
    Minimazer(std::function<double(double)> _func, double _a, double _b, double _eps=0.01, double _r=2.0) : func(_func) {
        a = _a; b = _b;
        eps = _eps;
        r= (_r > 1.0)? _r : 2.0;
    }
    dpair find_glob_min() {
        vector<pair<double,double> > vec;
        vec.push_back(dpair(a, func(a)));
        vec.push_back(dpair(b, func(b)));
        count = 2;
        double M=0;
        size_t k = 2;
        size_t t = 0;
        
        for (; ((vec[t+1].first) - (vec[t].first)) > eps;++k) {
            for (size_t i = 0; i < (k-1u); ++i) {
                double M_tmp = abs((vec[i+1].second-vec[i].second)/(vec[i+1].first - vec[i].first));
                 if (M_tmp > M) M = M_tmp;
            }

            double m = 1;
            if (M != 0) m = r * M;
            t = 0;
            double R= m * (vec[1].first - vec[0].first) + (pow((vec[1].second - vec[0].second), 2) / (m * (vec[1].first - vec[0].first)))-2*(vec[1].second+ vec[0].second);
            for (size_t i = 1; i < (k-1u); ++i) {
                double R_tmp = m*(vec[i+1].first - vec[i].first) + (pow((vec[i+1].second - vec[i].second),2) / (m*(vec[i+1].first - vec[i].first))) - 2 * (vec[i+1].second + vec[i].second);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            double x_t1 = (vec[t].first + vec[t + 1].first)/2 - (vec[t + 1].second - vec[t].second) / (2 * m);
            dpair t1_pair = dpair(x_t1, func(x_t1));
            ++count;
            vec.insert(std::lower_bound(vec.begin(), vec.end(), t1_pair, [](const dpair& a, const dpair& b) {return a.first <= b.first; }), t1_pair); //No need for sorting, only to insert
        }
        sol = vec[t]; solved = true;
        return vec[t];
    }
    
    void Show_info() {
        if (!solved) this->find_glob_min();
        std::cout << "Минимум функции = " << sol.second << " в точке x=" << sol.first << std::endl;
        std::cout << "Точность вычислений eps = "<< eps << std::endl;
        std::cout <<"Параметр метода r = "<< r << std::endl;
        std::cout << "Функция была посчитана " << count << " раз(а)" << std::endl;
    }

};

int main()
{
    setlocale(LC_ALL, "Russian");

    std::function<double(double)> Func1 = [](double x){ return sin(x); };

    Minimazer min1(Func1, 1, 4,0.01,2.0);
    min1.Show_info();
    std::cout << std::endl;

    Minimazer min2(Func1, 1, 4,0.01,1.1);
    min2.Show_info();
    std::cout << std::endl;

    return 0;
}