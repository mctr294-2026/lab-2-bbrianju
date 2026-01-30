#include "roots.hpp"
#include <cmath>
#include <limits>
#include <functional>
#include <utility>


bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    if (root == nullptr) return false; // No valid output location
    const double tol = 1e-6;
    const int max_iter = 1'000'000;
    
    // Evaluate function at endpoints to check for sign change
    double fa = f(a);
    double fb = f(b);
    
    if (fa * fb > 0) {
        // No root in [a, b]
        return false; 
    }

    // Iteratively shrink the interval [a, b] until convergence
    for (int i = 0; i < max_iter; ++i) {
        double mid = 0.5 * (a + b);
        double fm = f(mid);

        // Check convergence: function value close to zero or interval sufficiently small
        // or interval width within tolerance
        if (std::fabs(fm) < tol || std::fabs(b - a) <= tol) {
            *root = mid;
            return true;
        }
        
        // Determine which subinterval contains the root
        if (fa * fm < 0) {
            // Root is in [a, mid]
            b = mid;
            fb = fm;
        } else {
            // Root is in [mid, b]
            a = mid;
            fa = fm;
        }
    }

    // Maximum iterations reached; return midpoint as best estimate  
    *root = 0.5 * (a + b);
    return true;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    // Output pointer check
    if (root == nullptr) return false;
    
    const double tol = 1e-6;
    const int max_iter = 1'000'000;                     
    
    // Evaluate function at endpoints to check for sign change
    double fa = f(a);
    double fb = f(b);
    
    if (fa * fb > 0) return false; // No root in [a, b]
    if (fb - fa == 0) return false; // Prevent division by zero

    // Iteratively refine the estimate using the Regula Falsi formula
    for (int i = 0; i < max_iter; ++i) {
    
        // Compute the point where the secant line crosses the x-axis
        double c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        // Check convergence: function value close to zero or interval sufficiently small
        if (std::fabs(fc) < tol || std::fabs(b - a) <= tol) {
            *root = c;
            return true;    
        }

        // Update the interval [a, b] based on the sign of f(c)
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {            
            a = c;
            fa = fc;
        }
        
    }

    // Maximum iterations reached; return best estimate
    *root = (a * fb - b * fa) / (fb - fa);
    return true;
}
    
        

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {

    // Output pointer check                  
    if (root == nullptr) return false;

    // Parameters
    const double tol = 1e-6;
    const int max_iter = 1'000'000;

    if (a > b) std::swap(a, b); // Ensure a < b
    if (c < a || c > b) return false; // Initial guess not in [a, b]

    // Iteratively apply the Newton-Raphson formula
    for (int i = 0; i < max_iter; ++i) {
        double fc = f(c);

        // Check convergence
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

        // Fail if iteration leaves [a, b]
        if (!std::isfinite(c_new) || c_new < a || c_new > b) {
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

    // Normalize interval and validate initial guess
    if (a > b) std::swap(a, b);
    if (c < a || c > b) return false;


    double x0 = b;
    double x1 = c;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < max_iter; ++i) {
        // Check convergence
        if (std::fabs(f1) <= tol) {
            *root = x1;
            return true;
        }

        // Secant slope denominator 
        double denom = f1 - f0;
        if (std::fabs(denom) < 1e-12) {
            return false; // Prevent division by zero
        }   
    
        // Compute next approximation
        double x2 = x1 - f1 * (x1 - x0) / denom;

        // Fail if iteration leaves [a, b]
        if (!std::isfinite(x2) || x2 < a || x2 > b) {
            return false;
        }

        // Update for next iteration
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    return false; // Max iterations reached without convergence
}

