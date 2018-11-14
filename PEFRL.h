#ifndef PEFRL_H
#define PEFRL_H

#include <iostream>
#include "Body.h"
#include <vector>
using namespace std;

double xi = 0.1786178958448091;
double lambda = -0.2123418310626054;
double chi = -0.06626458266981849;

// gegeven de punten op hun beginposities, uitvoeren van de PEFRL methode
void PEFRL(vector<Body> bodies, double h, double tstart, double teind) {
	for (double j = tstart; j < teind; j += h) {
		for (unsigned int i = 0; i < bodies.size; ++i) {
			Vec r = bodies.at(i).getpos();
			Vec v = bodies.at(i).getvel();

			r += xi * h*v;
			bodies[i].setpos(r);

			v += 0.5*(1 - 2 * lambda)*h*calc_accel(bodies).at(i);

			r += chi * h*v;
			bodies[i].setpos(r);

			v += lambda * h*calc_accel(bodies).at(i);

			r += (1 - 2 * (chi + xi))*h*v;
			bodies[i].setpos(r);

			v += lambda * h*calc_accel(bodies).at(i);

			r += chi * h*v;
			bodies[i].setpos(r);

			v += 0.5*(1 - 2 * lambda)*h*calc_accel(bodies).at(i);

			r += xi * h*v;
			bodies[i].setpos(r);

			bodies[i].setvel(v);
		}
	}
}

#endif PEFRL_H
