#include "roots.hpp"
#include <cmath>
#include <limits>
#include <functional>


bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    if (root == nullptr) return false;
    const double tol = 1e-6;
    const int max_iter = 1'000'000;
    
    double fa = f(a);
    double fb = f(b);
    
    if (fa * fb > 0) return false; // No root in [a, b]

    for (int i = 0; i < max_iter; ++i) {
        double mid = 0.5 * (a + b);
        double fm = f(mid);

        if (std::fabs(fm) < tol || std::fabs(b - a) <= tol) {
            *root = mid;
            return true;
        }

        if (fa * fm < 0) {
            b = mid;
            fb = fm;
        } else {
            a = mid;
            fa = fm;
        }
    }
    *root = 0.5 * (a + b);
    return true;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    if (root == nullptr) return false;
    
    const double tol = 1e-6;
    const int max_iter = 1'000'000;                     
    
    double fa = f(a);
    double fb = f(b);
    
    if (fa * fb > 0) return false; // No root in [a, b]
    if (fb - fa == 0) return false;


    for (int i = 0; i < max_iter; ++i) {
    
    double c = (a * fb - b * fa) / (fb - fa);
    double fc = f(c);

        if (std::fabs(fc) < tol || std::fabs(b - a) <= tol) {
            *root = c;
            return true;    
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {            
            a = c;
            fa = fc;
        }
    }

    *root = (a * fb - b * fa) / (fb - fa);
    return true;
}
    
        

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    return false;
}

