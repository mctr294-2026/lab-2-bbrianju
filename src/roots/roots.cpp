#include "roots.hpp"
#include <cmath>
#include <limits>
#include <functional>
#include <utility>


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
    if (root == nullptr) return false;

    const double tol = 1e-6;
    const int max_iter = 1'000'000;

    if (a > b) std::swap(a, b);
    if (c < a || c > b) return false; // Initial guess not in [a, b]

    for (int i = 0; i < max_iter; ++i) {
        double fc = f(c);

        if (std::fabs(fc) <= tol) {
            *root = c;
            return true;
        }

        double gc = g(c);

        // Derivative must be usable
        if(!std::isfinite(gc) || std::fabs(gc) < 1e-12) {
            return false;
        }

        double c_new = c - fc / gc;

        // Bail if we get NaN/Inf
        if (!std::isfinite(c_new)) {
            return false;
        }

        c = c_new;

    }
        
    return false;

}



bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
     if (root == nullptr) return false;

    const double tol = 1e-6;
    const int max_iter = 1'000'000;

    double x0 = b;
    double x1 = c;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < max_iter; ++i) {
        if (std::fabs(f1) <= tol) {
            *root = x1;
            return true;
        }

        double denom = f1 - f0;
        if (std::fabs(denom) < 1e-12) {
            return false; // Prevent division by zero
        }   
    
        double x2 = x1 - f1 * (x1 - x0) / denom;

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    return false;
}

