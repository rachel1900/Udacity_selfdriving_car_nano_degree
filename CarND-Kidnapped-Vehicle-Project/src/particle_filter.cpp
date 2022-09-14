/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
//#include <Dense>
//#include <eigen3/Eigen/Dense>
#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
std::default_random_engine gen;  // creat a default_random_engine object gen
// using namespace Eigen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    int num_particles = 100;  // TODO: Set the number of particles
    
    // create a normal (Gaussian) distribution for noise of x, y and theta
    normal_distribution<double> noise_x_init(0, std[0]);
    normal_distribution<double> noise_y_init(0, std[1]);
    normal_distribution<double> noise_theta_init(0, std[2]);
    

    // create particle
    for (int i = 0; i < num_particles ; ++i) {
        Particle particle;  // creat a instance of structure Particale
        
        // Initialize all particles to first position (based on estimates of x, y, theta
        particle.id = i;
        particle.x = x;
        particle.y = y;
        particle.theta = theta;
        particle.weight=1.0;
        
        // add noise
        particle.x += noise_x_init(gen);
        particle.y += noise_y_init(gen);
        particle.theta += noise_theta_init(gen);
        
        // add particle to particles
        particles.push_back(particle);  // Set of current particles
        
    }
    // update signs of intialized
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
    //  (motion model)
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    
    // define normal distributions for sensor noise
    normal_distribution<double> noise_x(0, std_pos[0]);
    normal_distribution<double> noise_y(0, std_pos[1]);
    normal_distribution<double> noise_theta(0, std_pos[2]);
    
    // calculate new state
    for (int i = 0; i < num_particles ; ++i) {
        // calculate the new state (new position after time step)
        // check if there is rotation of the car about z-axis
        if (fabs(yaw_rate) < 0.00001) {   // there is no rotation about z-axis (2D movement)
              particles[i].x += velocity * delta_t * cos(particles[i].theta);
              particles[i].y += velocity * delta_t * sin(particles[i].theta);
            }
            else {   // there is rotation about z-axis (3D movement)
              particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
              particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
              particles[i].theta += yaw_rate * delta_t;
            }

            // add noise
            particles[i].x += noise_x(gen);
            particles[i].y += noise_y(gen);
            particles[i].theta += noise_theta(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
    
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    
    
    for (int i=0; i< observations.size(); i++){
        LandmarkObs landmark_obs = observations[i]; // read the Local (vehicle coords)  position of landmark observation
        
        // (initiallize distance value)
        // set the maxmium double as a start point to find the nearest neighbor
        double dist_min = std::numeric_limits<double>::max();
        //set -1 as initiallied map id
        int map_id = -1;
        
        for(int j =0 ; j < predicted.size(); j++){
            LandmarkObs landmark_predicted = predicted[i]; // predicted Vector of predicted landmark observations
            
            // calculate RSME of this position
            double current_dist = dist(landmark_obs.x, landmark_obs.y,
                                       landmark_predicted.x, landmark_predicted.y);
            
            // find the nearest neighbor
            if (current_dist < dist_min ){
                dist_min = current_dist; // update the nearest distance
                map_id = landmark_predicted.id; // update map_id
            }
        }
        
        // associat observation with a landmark identifier (nearest neighbor)
        observations[i].id=map_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

    for (int i=0; i<num_particles;i++){
        
        // Step1: translate VEHICLE'S coordinates into MAP'S coordinate
        // transform each observation marker from the vehicle's coordinates
        // to the map's coordinates, with respect to our particle.
        vector<LandmarkObs> transformed_obs;
        for(int j=0; j< observations.size();j++){
            // Homogenous Transformation Matrix
            double homogenous_tf_x = cos(particles[i].theta)* observations[j].x
                                    -sin(particles[i].theta)*observations[j].y
                                    +particles[i].x;
            double homogenous_tf_y = sin(particles[i].theta)* observations[j].x
                                    +cos(particles[i].theta)*observations[j].y
                                    +particles[i].y;
        
            transformed_obs.push_back(LandmarkObs{observations[j].id,homogenous_tf_x,homogenous_tf_y});
        }
 
    
        // Step 2: associate each transformed observation with a land mark identifier finde the nearest neighbor （only consider the landmarks within sensor_range）
        // for each particle create a vector to hold the map landmark that located within sensor range
        vector<LandmarkObs> landmark_selected;
        for (int m=0; m<map_landmarks.landmark_list.size(); m++){
            if (fabs(map_landmarks.landmark_list[m].x_f - particles[i].x) <= sensor_range
                &&
                fabs(map_landmarks.landmark_list[m].y_f - particles[i].y)){
                
                landmark_selected.push_back(LandmarkObs{map_landmarks.landmark_list[m].id_i,map_landmarks.landmark_list[m].x_f,
                    map_landmarks.landmark_list[m].y_f  });
            }
        }
        //identifier finde the nearest neighbor with help of function dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations)
        dataAssociation(landmark_selected,transformed_obs);
    
        // Step 3: update weight
        
        for (int k=0; k<transformed_obs.size();k++){
            for (int h=0; h<landmark_selected.size();h++){
                // check if the landmark_selected and transformed_obs are paired
                if (landmark_selected[h].id == transformed_obs[k].id){
                    // Multivariate-Gaussian probability
                    double sig_x = std_landmark[0]; //Landmark measurement uncertainty
                    double sig_y = std_landmark[1];
                    double gauss_norm;
                    double x_obs = transformed_obs[k].x;
                    double mu_x = landmark_selected[h].x;
                    double y_obs = transformed_obs[k].y;
                    double mu_y = landmark_selected[h].y;
            
                    gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
                    // calculate exponent
                    double exponent;
                    exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
                                   + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
                        
                    // calculate weight using normalization terms and exponent
                    double weight;
                    weight = gauss_norm * exp(-exponent);
                    // update weight
                    particles[i].weight *= weight;
                }
            }
        }
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    vector<Particle> resampled_particle;
    vector<double> weight_values;

    // get weihts of particles
    for (int i=0; i<particles.size();i++){
        weight_values.push_back(particles[i].weight);
    }
    
    double max_weight = *max_element(weight_values.begin(),weight_values.end());
    
    // Creating distributions.
    // uniform random distribution [0.0, 2.0*max_weight)
    // https://selfdriving5.github.io/udacity/Self-Driving%20Car%20Engineer%20v5.0.0%28us%29/Part%2002-Module%2001-Lesson%2004_Particle%20Filters/20.%20Resampling%20Wheel.html
    std::uniform_real_distribution<double>beta_list(0.0, 2.0*max_weight);
    
    // Generat random index for resample
    std::uniform_int_distribution<int>index_list(0, num_particles - 1);
    int index = index_list(gen);
    
    double beta = 0.0;
     // spin the resample wheel
     for (int i = 0; i < num_particles; i++) {
       beta += beta_list(gen);
       while (beta > weights[index]) {
         beta -= weights[index];
         index = (index + 1) % num_particles;
       }
       resampled_particle.push_back(particles[index]);
     }

     particles = resampled_particle;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
