#pragma once

#include "cgp/cgp.hpp"

// Particle structure 
struct particle_structure
{
	cgp::vec3 p0;  // Initial position
	float v0;  // Initial velocity
	float t0;      // Time of birth

	// Create a particle at its initial position
	particle_structure(cgp::vec3 p0, float v0, float creation_time=0);
	cgp::vec3 evaluate_position(float absolute_time) const;
};

struct particle_system_structure
{
	// Storage of the particles
	std::vector<particle_structure> particles; 

	// Create and add a new particle in the vector
	void create_new_particle(float current_time, float v0, cgp::vec3 p0);

	// Remove particles that are too old from the vector
	void remove_old_particles(float current_time);
};