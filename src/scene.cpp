#include "scene.hpp"


using namespace cgp;





void scene_structure::initialize()
{
	// Initialize the shapes of the scene
	// ***************************************** //
	N = 50;
	int Lfactor = 20;
	int Wfactor = 20;
	// Set the behavior of the camera and its initial position
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	//environment.camera.look_at({ 5.0f,-4.0f,2.0f }, { 0,0,0 });

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");

	// Create the ground plane
	ground.initialize(mesh_primitive_quadrangle({-2,-2,-1}, {2,-2,-1}, {2, 2,-1}, {-2, 2,-1}), "Ground");

	// send data to GPU and store it into a curve_drawable structure


	
	shape = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, Lfactor* N);
	for (int i = 0; i < N * N * Lfactor; i++){
		shape.position.at(i)=vec3(Wfactor * shape.position.at(i).x, Lfactor * shape.position.at(i).y, shape.position.at(i).z);
	}
	floor = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, Lfactor* N);
	for (int i = 0; i < N * N * Lfactor; i++){
		floor.position.at(i)=vec3(Wfactor * floor.position.at(i).x, Lfactor * floor.position.at(i).y,  atan_floor(floor.position.at(i).y)-floor_offset);//- 3.2);
	}
	environment.camera.look_at({ 5.0f,-4.0f,2.0f }, shape.position.at(N * N * Lfactor / 2 - N/2 ));
	//shape.position.at(N/2, 5*N/2, n/2)
	initial_position = shape.position;
	initial_time = timer.t;
	shape_visual.initialize(shape, "Deforming shape");
	shape_visual.shading.color = { 0.6f, 0.6f, 0.9f };
	floor_visual.initialize(floor, "Deforming floor");
	floor_visual.shading.color = { 0.6f, 0.6f, 0.0f };
		// Reset the color of the shape to white (only the texture image will be seen)
	shape_visual.shading.color = {1,1,1};

	// Load the image and associate the texture id to the structure
	//shape_visual.texture = opengl_load_texture_image("assets/squirrel.jpg");

}


void scene_structure::display()
{
	// Set the light to the current position of the camera
	environment.light = environment.camera.position();
	// Update the current elapsed time
	timer.update();


	// the general syntax to display a mesh is:
	//   draw(mesh_drawableName, environment);
	//     Note: scene is used to set the uniform parameters associated to the camera, light, etc. to the shader
	//draw(ground, environment);

	
	// conditional display of the global frame (set via the GUI)
	if (gui.display_frame)
		draw(global_frame, environment);


	draw(floor_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(floor_visual, environment, { 0,0,0 });

	draw(shape_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(shape_visual, environment, { 0,0,0 });
	
	evolve_shape();
	shape_visual.update_position(shape.position);

	// Recompute normals on the CPU (given the position and the connectivity currently in the mesh structure)
	shape.compute_normal();
	// Send updated normals on the GPU
	shape_visual.update_normal(shape.normal);

}


void scene_structure::display_gui()
{
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::SliderFloat("Time Scale", &timer.scale, 0.0f, 2.0f, "%.1f");
	ImGui::SliderFloat("Wind Strenght", &wind_str, 0.0f, 10.0f, "%.1f");
	ImGui::SliderFloat("Wind angle", &wind_angle, 0.0f, 3.14f, "%.01f");
	ImGui::SliderFloat("floor offset", &floor_offset, 0.0f, 10.f, "%.1f");
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	bool const restart = ImGui::Button("Restart");

	if (restart)
		initialize();
	if (restart)
		evolve_shape();
}

void scene_structure::evolve_shape()
{
    size_t const M = initial_position.size();

	//wind
	vec3 wind_dir = normalize(vec3(std::cos(wind_angle), std::sin(wind_angle), 0.0));
	vec3 wind_dir_2 = normalize(vec3(std::cos(wind_angle), -std::sin(wind_angle), 0.0));
	//wind_str = 3.f;

	float R = 0.007065f * std::pow(wind_str, 2.5) / 2; //0.06;
	float w = 9.81f * std::sqrt(2.f/3.f) / wind_str ; //4.0;

	float K = 9.81f * 2 / 3 / wind_str / wind_str;   //10.0;
	vec3 K_vec = K * wind_dir;
	float lambda = 1.5;

	float K_y = 1.0;
	float K_z = 1.0;
	//std::cout << "R = " << R << " w = " << w << " K = " << K <<"\n";
    for(size_t k=0; k<M; ++k)
    {
		vec3 const& p0 = initial_position[k];
        vec3& p        = shape.position[k];
		float gamma = 0;
		if (k>N && k < M-N) gamma = (floor.position[k+N].z - floor.position[k-N].z)/(floor.position[k+N].y - floor.position[k-N].y);
		else if (k > M-N-1) gamma = (floor.position[M-1].z - floor.position[M-1-N].z)/(floor.position[M-1].y - floor.position[M-1-N].y);
		float depth =  -1*floor.position[k].z;
		float K_depth = K / std::sqrt( std::tanh(K * depth));
		float w_depth = 0;

		float phi = - w * timer.t + K_depth * dot(wind_dir, p0) ;//K_integration(K ,dot(wind_dir, p0), &valley_floor);
		//float alpha = -1 * K_depth * dot(wind_dir_2, p0);
		float alpha = 0.23 * std::exp(-0.1 * K_depth * depth ) ;
		float Sx = 1 / (1 - std::exp(-1 * K_y * depth));
		float Sz = Sx * (1 - std::exp(-1 * K_z * depth));

		float alpha2 = std::asin(sin(gamma)*std::exp(-K_depth*depth));
		
        //p.x = p0.x + R * dot(wind_dir, vec3(1,0,0)) * std::sin(K * p0.y - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));
		//p.y = p0.y + R  * std::sin( K_depth * dot(wind_dir, p0) - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));
		//p.z = p0.z - R * std::cos(K_depth * dot(wind_dir, p0)  - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));


		//p.y= p0.y + 2*R * std::cos(phi + alpha - lambda * (p.z - p0.z) * (timer.t - initial_time));
		//p.z= p0.z +3*R * std::sin(phi + alpha - lambda * (p.z - p0.z) * (timer.t - initial_time));

		p.y = p0.y + R * std::cos(alpha2) * Sx * std::sin(phi) +
			 std::sin(alpha2) * Sz * std::cos(phi);
		p.z = p0.z + R *  std::cos(alpha2) * Sz * std::sin(phi) +
			 std::sin(alpha2) * Sx * std::cos(phi);
			
    }
}

float K_integration(float K ,float x0, float (*h)(float, float, float)){
	float sum = 0;
	
	for (int i = 0; i < 100; i++){
		sum+= K / std::sqrt(K * h(i * x0/100, 2.0, 0.23)) * x0 /100;
	}
	return sum;
}

float sloped_floor(float x, float limit, float slope){
	if (x < limit){
		return slope * limit;
	}
	else {
		return  slope * (x);
	}
}

float atan_floor(float x){
	return 0.5*(std::atan(40*(x-0.4))-1.7);
}

float valley_floor(float x, float position, float width){
	
	if ( x >  position - 3 * width) if ( x < position + 3 * width) return 3.0;
	return 0.6;
}

float valley_wall(float x, float position, float width){
	
	if ( x >  position ) return 13.0;
	return 0.6;
}

