#include "particles.hpp"

using namespace cgp;

particle_structure::particle_structure(vec3 p, float v, float creation_time)
{
	p0 = p;
	v0 = v;
	t0 = creation_time;
}

vec3 particle_structure::evaluate_position(float absolute_time) const
{
	//vec3 const g = { 0,0,-9.81f };      // gravity constant
	float const t = absolute_time - t0; // local time elapsed since the particle creation

	// TO DO: Modify this computation to model the bouncing effect
	// **************************************************************** //
	vec3 p;
	p.z = p0.z - 0.5f * 9.81 * t * t *0.1 ;
	p.y = p0.y + v0*5.0*t;
	p.x = p0.x;

	return p;
}

void particle_system_structure::create_new_particle(float t0, float v0, vec3 p0)
{
	particles.push_back(particle_structure(p0, v0, t0));
}

void particle_system_structure::remove_old_particles(float t)
{
	float const max_time = 1.0f;

	// Loop over all active particles
	for (auto it = particles.begin(); it != particles.end();)
	{
		// if a particle is too old, remove it
		if (t - it->t0 > max_time) {

			// remove particle
			it = particles.erase(it);

		}

		// Go to the next particle if we are not already on the last one
		if (it != particles.end()) {
			++it;
		}
	}
}