// C++ program to determine the minimum length of the supersonic nozzle
// based on Method of Characteristics
// Generates x-y plot and the nozzle x-y coordinates is saved in a file "nozzle_xy.txt"

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>

#include "matplotlibcpp.h"                  // python plot library
namespace plt = matplotlibcpp;

double PM_function(double &, double &);
double mach_angle_fn(double &, double &);
double inverse_PM_function(double &, double &);

int main() {
    std::cout << std::fixed << std::setprecision(4);        

    /*-----------------------------USER INPUT-----------------------------*/
    double gamma {1.4};                // Specify heat ratio of air
    double radius {35};                // Specify half-throat radius "cm"
    double mach_exit {2.0};            // Specify exit Mach number
    unsigned int num_lines {25};       // Specify number of left running characteristic lines (atleast 2) == num of points for nozzle wall contour + 1
    /*--------------------------------------------------------------------*/

    if (num_lines < 2) {
        std::cerr << "Insufficient Data. Number of lines should be greater than or equal to 2." << std::endl;
        std::cin.ignore();
    }

    else {
            double x_org{ 0 }, y_org{ radius };     // position of half-throat section
            double nu{}, mu{}, theta{}, theta_wall_max{}, Ma; // nu - prandtl-meyer angle; mu - mach angle
            double K_left{}, K_right{}, slope_left{}, slope_right{}, rhs_l{}, rhs_r{};

            theta_wall_max = PM_function(mach_exit, gamma) / 2.0;  // maximum angle at the throat

            std::vector<double> x_0, y_0, theta_0, nu_0, Ma_0, mu_0;
            std::vector<double> slope_right_0, slope_left_0, K_right_0, K_left_0;

            std::vector<double> x_wall, y_wall;
            std::vector<double> x_s(2), y_s(2), x_w(2), y_w(2);     // for plotting "moc" with python matplotlib library

            for (int i{ 0 }; i < num_lines; i++) {
                x_0.push_back(x_org);
                y_0.push_back(y_org);         // y-coordinate is at 1

                x_s.at(0) = x_org;
                y_s.at(0) = y_org;

                if (i == 0)
                    theta = theta_wall_max / num_lines;
                else
                    theta = (theta_wall_max / num_lines) + theta_0[i - 1];

                theta_0.push_back(theta);
                nu = theta;
                nu_0.push_back(nu);
                Ma = inverse_PM_function(nu, gamma);
                mu = mach_angle_fn(Ma, gamma);
                Ma_0.push_back(Ma);
                mu_0.push_back(mu);
                K_left = theta + nu;
                K_right = theta - nu;
                K_left_0.push_back(K_left);
                K_right_0.push_back(K_right);
            }

            int n = num_lines;
            int r = n + 1;
            int s = r;

            std::vector<double> x_p, y_p, theta_p, nu_p, Ma_p, mu_p;
            std::vector<double> slope_right_p, slope_left_p, K_right_p, K_left_p;
            std::vector<double> x_n(s), y_n(s), theta_n(s), nu_n(s), Ma_n(s), mu_n(s), K_left_n(s), K_right_n(s);
            double x{ 0.0 }, y{ 0.0 };

            for (int i{ 0 }; i < n; i++) {
                for (int j{ 0 }; j < r; j++) {

                    if (j == 0) {
                        theta = 0;      // line of symmetry condition
                        y = 0;          // line of symmetry condition

                        if (i == 0)
                            K_left = K_left_0[i];

                        else
                            K_left = K_left_n[j + 1];

                        nu = K_left - theta;
                        K_right = theta - nu;
                        Ma = inverse_PM_function(nu, gamma);
                        mu = mach_angle_fn(Ma, gamma);

                        if (i == 0)
                            rhs_l = ((theta + theta_0[i]) - (mu + mu_0[i])) * (M_PI / 180) / 2.0;

                        else
                            rhs_l = ((theta + theta_n[j + 1]) - (mu + mu_n[j + 1])) * (M_PI / 180) / 2.0;

                        slope_left = tan(rhs_l);
                        slope_right = 0;

                        if (i == 0) {
                            x = x_0[i + 1] + ((y - y_0[i + 1]) / slope_left);     // Finding the coordinates of a point on line of symmetry
                        }
                        else {
                            x = x_n[j + 1] + ((y - y_n[j + 1]) / slope_left);     // Finding the coordinates of a point on line of symmetry
                        }

                        x_n[j] = x;
                        y_n[j] = y;
                        theta_n[j] = theta;
                        nu_n[j] = nu;
                        Ma_n[j] = Ma;
                        mu_n[j] = mu;
                        K_left_n[j] = K_left;
                        K_right_n[j] = K_right;

                        x_s.at(1) = x;
                        y_s.at(1) = y;
                        x_w.at(0) = x;
                        y_w.at(0) = y;

                        plt::plot(x_s, y_s, "-g");      // plot x and y of centerline section
                    }

                    else if (j == r - 1) {              // last ch. point -- at the wall
                        theta = theta_p[j - 1];         // wave cancellation condition
                        nu = K_left - theta;
                        K_right = theta - nu;
                        Ma = inverse_PM_function(nu, gamma);
                        mu = mach_angle_fn(Ma, gamma);
                        double E{};

                        if (i == 0) {
                            rhs_r = ((theta + mu_n[j - 1])) * (M_PI / 180);
                            slope_right = tan(rhs_r);

                            rhs_l = (theta + theta_wall_max) * (M_PI / 180) / 2.0;
                            slope_left = tan(rhs_l);

                            E = slope_right / slope_left;

                            y = (1 / (1 - (E))) * (y_n[j - 1] - (y_0[i] * (E))
                                + (slope_right * (x_0[i] - x_n[j - 1])));

                            x = x_n[j - 1] + ((y - y_n[j - 1]) / slope_right);

                        }
                        else {

                            rhs_r = ((theta + mu)) * (M_PI / 180);
                            slope_right = tan(rhs_r);

                            rhs_l = (theta + theta_n[j + 1]) * (M_PI / 180) / 2.0;
                            slope_left = tan(rhs_l);

                            E = slope_right / slope_left;

                            y = (1 / (1 - (E))) * (y_n[j - 1] - (y_n[j + 1] * (E))
                                + (slope_right * (x_n[j + 1] - x_n[j - 1])));

                            x = x_n[j - 1] + ((y - y_n[j - 1]) / slope_right);

                        }

                        x_n[j] = x;
                        y_n[j] = y;
                        theta_n[j] = theta;
                        nu_n[j] = nu;
                        Ma_n[j] = Ma;
                        mu_n[j] = mu;
                        K_left_n[j] = K_left;
                        K_right_n[j] = K_right;

                        x_wall.push_back(x);
                        y_wall.push_back(y);

                        x_w.at(1) = x;
                        y_w.at(1) = y;

                        plt::plot(x_w, y_w, "-b");      // plot x and y of wall section

                    }
                    else {
                        if (i == 0) {
                            K_left = K_left_0[j];
                            K_right = K_right_p[j - 1];
                        }
                        else {
                            K_left = K_left_n[j + 1];
                            K_right = K_right_n[j - 1];
                        }

                        theta = (K_left + K_right) / 2;
                        nu = K_left - theta;
                        Ma = inverse_PM_function(nu, gamma);
                        mu = mach_angle_fn(Ma, gamma);

                        if (i == 0) {
                            rhs_l = ((theta + theta_0[j]) - (mu + mu_0[j])) * (M_PI / 180) / 2.0;
                            rhs_r = ((theta + theta_p[j - 1]) + (mu + mu_p[j - 1])) * (M_PI / 180) / 2.0;
                        }

                        else {
                            rhs_l = ((theta + theta_n[j + 1]) - (mu + mu_n[j + 1])) * (M_PI / 180) / 2.0;
                            rhs_r = ((theta + theta_n[j - 1]) + (mu + mu_n[j - 1])) * (M_PI / 180) / 2.0;
                        }

                        slope_right = tan(rhs_r);
                        slope_left = tan(rhs_l);

                        double A{}, B{}, C{}, D{};

                        if (i == 0) {
                            A = (1 / (1 - (slope_right / slope_left)));
                            B = y_p[j - 1] - (slope_right * x_p[j - 1]);
                            C = (y_0[i + 1] * slope_right / slope_left);
                            D = x_0[i + 1] * slope_right;
                            y = A * (B - C + D);
                            x = x_p[j - 1] + ((y - y_p[j - 1]) / slope_right);
                        }
                        else {
                            A = (1 / (1 - (slope_right / slope_left)));
                            B = y_n[j - 1] - (slope_right * x_n[j - 1]);
                            C = (y_n[j + 1] * slope_right / slope_left);
                            D = x_n[j + 1] * slope_right;
                            y = A * (B - C + D);
                            x = x_n[j - 1] + ((y - y_n[j - 1]) / slope_right);
                        }

                        x_n[j] = x;
                        y_n[j] = y;
                        theta_n[j] = theta;
                        nu_n[j] = nu;
                        Ma_n[j] = Ma;
                        mu_n[j] = mu;
                        K_left_n[j] = K_left;
                        K_right_n[j] = K_right;

                    }

                    x_p.push_back(x);
                    y_p.push_back(y);
                    theta_p.push_back(theta);
                    nu_p.push_back(nu);
                    Ma_p.push_back(Ma);
                    mu_p.push_back(mu);
                    K_left_p.push_back(K_left);
                    K_right_p.push_back(K_right);
                    slope_left_p.push_back(slope_left);
                    slope_right_p.push_back(slope_right);

                }
                r--;
                s--;

                if (r < 2)
                    break;
            }

            // Write the nozzle xy coordinates to a ".txt" file
            std::ofstream ofile;
            ofile.open("nozzle_xy.txt");
            ofile << std::fixed << std::setprecision(4);

            ofile << std::left << std::setw(10) << "x\t"
                << std::setw(10) << "y\t"
                << std::endl;

            ofile << std::left << std::setw(12) << x_org
                << std::setw(10) << y_org
                << std::endl;

            for (size_t i{ 0 }; i < x_wall.size(); i++) {
                ofile << std::left << std::setw(10) << x_wall[i] << '\t'
                    << std::setw(10) << y_wall[i] << '\t'
                    << std::endl;
            }
            ofile.close();

            // Plot properties
            double x_max = x_w.at(1) + 5;
            double y_max = y_w.at(1) + 5;
            plt::title("Method of Characteristics");
            plt::xlabel("Centerline (cm)");
            plt::ylabel("Radius (cm)");
            plt::xlim(0.0, x_max);
            plt::ylim(0.0, y_max);
            plt::save("moc.jpg");           // plot file name
            plt::show();
    }
    
    return 0;
}

// Calculate Prandtl-Meyer Angle from Mach Number
double PM_function(double &mach, double &gamma) {
    double PM_angle{};
    PM_angle = ((sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt(((pow(mach, 2) - 1) * (gamma - 1) / (gamma + 1)))))
        - (atan(sqrt(pow(mach, 2) - 1)))) * 180 / M_PI;
    return PM_angle;    // in degrees
}

// Calculate Mach Angle for a corresponding Mach Number
double mach_angle_fn(double &mach, double &gamma) {
    double mach_ang{};
    mach_ang = (asin(1.0 / mach)) * 180 / M_PI;
    return mach_ang;    // in degrees
}

// Calculate Mach Number from Prandtl-Meyer Angle
// Based on "Inversion of the Prandtl-Meyer Relation" by I.M. Hall (Aeronautical Journal, September 1975)
double inverse_PM_function(double &theta, double &gamma) {
    double beta{}, F_beta{}, beta_n{ 0.0 };
    double mach_fi;         // final mach value for a corresponding PM Angle
    double mach{ 1.001 };   // initial guess
    double diff{ 1 }, lambda{}, nu_m{};
    int c{ 1 };

    lambda = sqrt((gamma - 1) / (gamma + 1));
    nu_m = theta * M_PI / 180;

    while (c != 0) {
        beta = sqrt(pow(mach, 2) - 1);
        F_beta = ((1.0 / lambda) * tan(lambda * (nu_m + atan(beta)))) - beta;
        beta_n = beta + (((1.0 + pow(beta, 2.0)) * F_beta) / ((1.0 - pow(lambda, 2.0)) * pow(beta, 2.0)));
        mach_fi = static_cast<float>(sqrt(pow(beta_n, 2) + 1.0));
        diff = mach - mach_fi;
        mach = mach_fi;

        if (diff == 0)
            c = 0;
    }
    return mach_fi;
}