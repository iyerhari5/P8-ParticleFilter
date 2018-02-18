/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;

	// This line creates normal (Gaussian) distributions for x,y, and theta.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	default_random_engine gen;

	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);

		Particle p;
		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = 1;

		particles.push_back(p);
	}

	is_initialized = true;
	//cout << particles[0].x << ", " << particles[0].y << ", " << particles[0].weight << endl;
	cout << "Initialized!!" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	double xf, yf, thetaf;
	default_random_engine gen;

	float thetadot_deltat = yaw_rate*delta_t;
	for (int i = 0; i < num_particles; ++i) {
				
		if (fabs(yaw_rate) < 0.00001) {
			xf = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			yf = particles[i].y + velocity * delta_t * sin(particles[i].theta);
		}
		else
		{
			xf = particles[i].x + velocity / yaw_rate*(sin(particles[i].theta + thetadot_deltat) - sin(particles[i].theta));
			yf = particles[i].y + velocity / yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + thetadot_deltat));
			particles[i].theta = particles[i].theta + thetadot_deltat;
		}

		//add sensor noise to the positions
		normal_distribution<double> dist_x(xf, std_pos[0]);
		particles[i].x = dist_x(gen);

		normal_distribution<double> dist_y(yf, std_pos[1]);
		particles[i].y = dist_y(gen);
	}

	//cout << "Predicted!!" << endl;
	//cout << particles[0].x << ", " << particles[0].y << ", "<< particles[0].weight  << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	
	for (int i = 0; i < observations.size(); i++) {
		double minDist = 100000000;
		int minDistId = 100000;
		double distance;

		for (int j = 0; j < predicted.size(); j++) {

			distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if (distance < minDist) {
				minDist = distance;
				minDistId = j;
			}
		}
		//assign the closest object
		observations[i].id = minDistId;// predicted[minDistId].id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	//cout << " Inside update weights" << endl;
	std::vector<LandmarkObs> predicted;
	float xP, yP;
	double distance;
	float theta;

	double sumWeights = 0;

	std::vector<LandmarkObs> observationsMap;
	std::vector<int> selectedLandmarks;
	for (int p = 0; p < num_particles; p++) {

		//convert the observations to the map coordinates
		observationsMap.clear();
		theta = particles[p].theta;
		for (int i = 0; i < observations.size(); i++) {
			xP = particles[p].x + cos(theta)*observations[i].x - sin(theta)*observations[i].y;
			yP = particles[p].y + sin(theta)*observations[i].x + cos(theta)*observations[i].y;
			LandmarkObs landmark;
			landmark.x = xP;
			landmark.y = yP;
			landmark.id = observations[i].id;
			observationsMap.push_back(landmark);
		}


		predicted.clear();
		for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
			xP = map_landmarks.landmark_list[i].x_f;
			yP = map_landmarks.landmark_list[i].y_f;

			distance = dist(xP, yP, particles[p].x, particles[p].y);
			if (distance < sensor_range)
			{
				LandmarkObs landmark;
				landmark.x = xP;
				landmark.y = yP;
				landmark.id = map_landmarks.landmark_list[i].id_i;
				predicted.push_back(landmark);
			}
		}

		//Now we have all the predicted landmark points within sensor range of this particle
		//We will find the data associations now
		dataAssociation(predicted, observationsMap);
		/*
		particles[p].associations.clear();
		for (int m = 0; m < observationsMap.size(); m++) {
			particles[p].associations.push_back(predicted[observationsMap[m].id].id);
		}*/

		//Update the weight for this particle
		double two_std_x_sqr = 2 * std_landmark[0] * std_landmark[0];
		double two_std_y_sqr = 2 * std_landmark[1] * std_landmark[1];
		double denom = 1.0 / (2 * M_PI*std_landmark[0] * std_landmark[1]);

		particles[p].weight = 1;
		for (int j = 0; j < observationsMap.size(); j++) {

			double nearestX = predicted[observationsMap[j].id].x;
			double nearestY = predicted[observationsMap[j].id].y;
			double exponent = (observationsMap[j].x - nearestX)*(observationsMap[j].x - nearestX) / two_std_x_sqr +
				(observationsMap[j].y - nearestY)*(observationsMap[j].y - nearestY) / two_std_y_sqr;
			double prob = denom*exp(-exponent);
			particles[p].weight *= prob;
		}//for each observation

		//add the weight to the sum
		sumWeights += particles[p].weight;

	}//for each particle

	//Now normalize the weights to be between 0 and 1
	for (int p = 0; p < num_particles; p++)
		particles[p].weight /= sumWeights;

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<> dn(0,num_particles-1);
	int index = dn(gen);
	double beta = 0;

	//find the maximum weight among the particles
	float maxWeight = 0;
	for (int i = 0; i < num_particles; i++)
		if (particles[i].weight > maxWeight)
			maxWeight = particles[i].weight;
	
	std::vector<Particle> new_particles;
	uniform_real_distribution<double> ud(0, maxWeight);
	for (int i = 0; i < num_particles; i++) {
		beta += ud(gen);
		while (beta > particles[index].weight) {
			beta -= particles[index].weight;
			index = (index + 1) % (num_particles - 1);
		}
		new_particles.push_back(particles[index]);
	}

	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
