#include "scene.hpp"
#include <random>


using namespace cgp;

void scene_structure::initialize()
{
	// Initialize the shapes of the scene
	// ***************************************** //
	//int N = 100;
	
	int Lfactor = 20;
	int Wfactor = 50;

	y_collision.resize(N);
	y_speed_waterline.resize(N);

	int N2 = Lfactor * N;
	// Set the behavior of the camera and its initial position
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	//environment.camera.look_at({ 5.0f,-4.0f,2.0f }, { 0,0,0 });

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");

	// Create the ground plane
	ground.initialize(mesh_primitive_quadrangle({-2,-2,-1}, {2,-2,-1}, {2, 2,-1}, {-2, 2,-1}), "Ground");

	// send data to GPU and store it into a curve_drawable structure

	shape = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, N2);
	shape_waterline = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, N);
	floor = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, N2);

	for (int i = 0; i < N * N2; i++){ // adjusting floor position
		floor.position.at(i)=vec3(Wfactor * floor.position.at(i).x, Lfactor * floor.position.at(i).y,
								 (sloped_floor_xy(Lfactor * floor.position.at(i).y, Wfactor * floor.position.at(i).x, 0.23) -floor_offset)* noise_perlin(floor.position.at(i).x, octave, persistance, gain) );
	}

	// detecting waterline
	float epsilon = 0.05;
	int j,k;

	for (int i = 0; i < N; i++) y_collision[i] = Lfactor; // initializing waterline
	for (int i = 0; i < N; i++) y_speed_waterline[i] = 0.0; // initializing waterline speed

	for (int i = 0; i < N * N2; i++){ // getting waterline value
		k = i / N;
		j = i % N;

		if (std::abs(floor.position.at(i).z) < epsilon && floor.position.at(i).y < y_collision[j]) {
			y_collision[j] = floor.position.at(i).y;
		}
	}
	
	// adjusting shape position
	for (int i = 0; i < N * N2; i++) shape.position.at(i)=vec3(Wfactor * shape.position.at(i).x, shape.position.at(i).y * y_collision[j] , shape.position.at(i).z);
	// adjusting waterline position
	for (int i = 0; i < N * N; i++){
		shape_waterline.position.at(i) = vec3(Wfactor * shape_waterline.position.at(i).x, 
		shape_waterline.position.at(i).y*6 + y_collision[j] - 4.5 ,
		shape_waterline.position.at(i).z);
	}

	environment.camera.look_at({ 5.0f,-4.0f,2.0f }, shape.position.at(N * N * Lfactor / 2 - N/2 ));
	initial_position = shape.position;
	initial_waterline = shape_waterline.position;
	//for (int l = 0; l < N; l++) y_collision[l] = initial_position[initial_position.size()-1].y;
	initial_time = timer.t;
	shape_visual.initialize(shape, "Deforming shape");
	shape_visual.shading.color = { 0.6f, 0.6f, 0.9f };
	shape_waterline_visual.initialize(shape_waterline, "Deforming waterline");
	shape_waterline_visual.shading.color = { 0.6f, 0.6f, 0.9f };
	floor_visual.initialize(floor, "Deforming floor");
	floor_visual.shading.color = { 0.6f, 0.6f, 0.0f };

	float sphere_radius = 0.05f;
	sphere.initialize(mesh_primitive_sphere(sphere_radius), "Sphere");
	sphere.shading.color = { 1.0f, 1.0f, 1.0f };

	// Reset the color of the shape to white (only the texture image will be seen)
	shape_visual.shading.color = {1,1,1};
	shape_waterline_visual.shading.color = {1,1,1};
	// Load the image and associate the texture id to the structure
	if (texturesOn){
		shape_visual.texture = opengl_load_texture_image("assets/sea.jpg");
		shape_waterline_visual.texture = opengl_load_texture_image("assets/sea.jpg");
		floor_visual.texture = opengl_load_texture_image("assets/sand.jpg");
	}

}


void scene_structure::display()
{
	// Set the light to the current position of the camera
	environment.light = environment.camera.position();
	// Update the current elapsed time
	timer.update();

	// conditional display of the global frame (set via the GUI)
	if (gui.display_frame)
		draw(global_frame, environment);


	draw(floor_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(floor_visual, environment, { 0,0,0 });

	draw(shape_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(shape_visual, environment, { 0,0,0 });

	draw(shape_waterline_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(shape_waterline_visual, environment, { 0,0,0 });
	
	if (foamOn) particle_system.remove_old_particles(timer.t);

	// Evaluate the positions and display the particles
	int const N = particle_system.particles.size();
	for (int k = 0; k < N; ++k)
	{
		// Current particle
		particle_structure& particle = particle_system.particles[k];

		// Evaluate the current position of the particle
		vec3 const p = particle.evaluate_position(timer.t);

		// Display the particle as a sphere
		sphere.transform.translation = p;
		draw(sphere, environment);
	}

	evolve_shape();
	if (foamOn) evolve_foam(timer.t);
	shape_visual.update_position(shape.position);
	shape_waterline_visual.update_position(shape_waterline.position);
	// Recompute normals on the CPU (given the position and the connectivity currently in the mesh structure)
	shape.compute_normal();
	shape_waterline.compute_normal();
	// Send updated normals on the GPU
	shape_visual.update_normal(shape.normal);
	shape_waterline_visual.update_normal(shape_waterline.normal);

}


void scene_structure::display_gui()
{
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Enable kludge", &kludge);
	ImGui::Checkbox("Enable textures", &texturesOn);
	ImGui::Checkbox("Enable foam (slow)", &foamOn);
	ImGui::SliderFloat("Time Scale", &timer.scale, 0.0f, 2.0f, "%.1f");
	ImGui::SliderFloat("Wind Strenght", &wind_str, 0.1f, 10.0f, "%.1f");
	ImGui::SliderFloat("Wind angle", &wind_angle, 0.7f, 1.7f, "%.01f");
	ImGui::SliderFloat("floor offset", &floor_offset, 0.0f, 10.f, "%.1f");
	ImGui::SliderFloat("Wave height", &K_var, 0.0f, 8.f, "%.1f");
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	bool const restart = ImGui::Button("Restart");

	if (restart)
		initialize();
	if (restart)
		evolve_shape();
}

void scene_structure::evolve_shape()
{
    int M = initial_position.size();

	//int N = 100;
	//wind
	vec3 wind_dir = normalize(vec3(std::cos(wind_angle), std::sin(wind_angle), 0.0));
	vec3 wind_dir_2;
	if (wind_angle < 1.6 ) wind_dir_2 = normalize(vec3(std::cos(wind_angle), std::sin(wind_angle), 0.0));
	else wind_dir_2 = normalize(vec3( std::cos(wind_angle), -1*std::sin(wind_angle), 0.0));

	float R = 0.007065f * std::pow(wind_str, 2.5) / 2; //0.06;
	float w = 9.81f * std::sqrt(2.f/3.f) / wind_str ; //4.0;
	float K = 9.81f * 2 / 3 / wind_str / wind_str;   //10.0;

	float lambda = 2.5;

	float K_o;
	float K_y = 1.0;
	float K_z = 20.0;
	float epsilon =  0.05;

	int i,j;

    for(size_t k=0; k<M; ++k)
    {
       
		i = k / N;
		j = k % N;

		vec3 const& p0 = initial_position[k];
        vec3& p        = shape.position[k];
		float depth =  std::sqrt(floor.position[k].z * floor.position[k].z  + 0.1); // to avoid zero values for depth
		float K_depth = K / std::sqrt( std::tanh(K * depth));

		float phi = - w  * timer.t + K_depth * dot(wind_dir, p0);
		if (kludge) phi -= lambda * (p.z - p0.z) * ((timer.t - initial_time) - static_cast<int>(timer.t - initial_time)%static_cast<int>(5*timer.t) );
		
		float phi2 = - w  * timer.t + K_depth * dot(wind_dir_2, p0);;
		float gamma = 0.0;
		if (k > N) gamma = (floor.position[k].z - floor.position[k - N].z)/(floor.position[k].y - floor.position[k - N].y);
		 
		float alpha = std::asin(std::sin(gamma) * std::exp(-0.1 * K_depth * depth )) ;
		float alpha2 = std::asin(std::sin(gamma) * std::exp(-0.1 * K_depth * depth ));

		float Sx = 1 / (1 - std::exp(-1 * K_y * depth));
		float Sz = Sx * (1 - std::exp(-1 * K_z * depth));

		p.y = p0.y + R * std::cos(alpha) * Sx * std::sin(phi) + std::sin(alpha) * Sz * std::cos(phi);
		p.z = p0.z + R *  std::cos(alpha) * Sz * std::sin(phi) + std::sin(alpha) * Sx * std::cos(phi);

		float train_profile = K_var * GaussianGate(dot(cross(wind_dir, vec3(0,0,1)), p0), 2.5, 20.0, 35.0)  * ( 1.0 + 0.3 * noise_perlin(p0.x,1,0.3,1.0));
		p.z *= train_profile ;

		// deforming the waterline 
		if (i ==  M/N - 1 - 3 * N){ 
			for(int l = 0; l < N; l++){
				shape_waterline.position[l * N + j].y = initial_waterline[j].y + 
														2 * K_var * (shape.position[l * M/N + j].z * (1 + 0.001 / train_profile) - initial_position[l * M/N + j].z) +
														2 * K_var * (0.1 + std::abs(wind_angle - 1.6)) * (R *  std::cos(alpha) * Sz * std::sin(phi2) + std::sin(alpha) * Sx * std::cos(phi2))+
														0.5 * noise_perlin(p0.x,1,0.3,2.0);

				shape_waterline.position[l * N + j].z = 0.05  + (sloped_floor_xy(shape_waterline.position[l * N + j].y, shape_waterline.position[l * N + j].x, 0.23) - floor_offset) ;//*  noise_perlin(shape_waterline.position[l * N + j].x, octave, persistance, gain);
			}
		}
	}
		
}

void scene_structure::create_foam_train(float t0, int k, float foam_th){
	float d = 0;
	if (k>N) d = (shape.position[k].z - shape.position[k-N].z)/(shape.position[k].y - shape.position[k-N].y);
	//std::cout<<d<<"\n";
	std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(0, 1);
	if (d > foam_th && distr(eng) < 0.001){
		particle_system.create_new_particle(t0, d, shape.position[k]);
		// for (int i = 0; i < 50; i++){
		// 	vec3 p0 = shape.position[k];
		// 	p0.x+=0.01*i;
		// 	particle_system.create_new_particle(t0, d, p0);
		// }
	}
}

void scene_structure::evolve_foam(float t0){
	size_t const M = initial_position.size();
	float foam_th = 0.15;
	//std::cout<<particle_system.particles.size()<<"BEFORE\n";
	int nb_train = 0;
	for(size_t k=0; k<M; ++k){
		create_foam_train(t0, k, foam_th);
	}
	//std::cout<<particle_system.particles.size()<<"AFTER\n";
}

float K_integration(float K ,float x0, float (*h)(float, float, float)){
	float sum = 0;
	
	for (int i = 0; i < 100; i++){
		sum+= K / std::sqrt(K * h(i * x0/100, 2.0, 0.23)) * x0 /100;
	}
	return sum;
}

float sloped_floor(float x, float limit, float slope){
	/* if (x < limit){
		return slope * limit;
	}
	else { */
		return  slope * (x);
	//}
}

float sloped_floor_xy(float x, float y, float slope){
	return slope * (x - std::sin(y/5));
}

float valley_floor(float x, float position, float width){
	
	if ( x >  position - 3 * width) if ( x < position + 3 * width) return 3.0;
	return 0.6;
}

float valley_wall(float x, float position, float width){
	
	if ( x >  4* position ) return 4.0;
	return 0.6;
}

float atan_floor(float x, float dist_from_shore, float steepness ){
	return 0.5*(std::atan(steepness*(x-dist_from_shore))-1.7);
	//return 0.3 * (std::atan(20 * (x - 10)) - 1.6) + 0.07*x;
}

float Gaussian(float x, float std, float avg ){
	return std::exp(-1*(x-avg)*(x-avg)/2/std/std) / std::sqrt(2*M_1_PI) / std;
}

float GaussianGate(float x, float std, float start, float end ){
	int n = static_cast<int>((end - start) / std / 2) + 1;
	float gate = 0.0;
	for (int k = 0; k < n; k++) gate += Gaussian(x, std, start + 2 * k * std);
	return gate / n ;
}



