#include "Simulator.h"
#include "File.h"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <iomanip>

void Simulator::Reset(Simulation *_simulation, SimulatorDesc *const _simulator_desc)
{
	simulation = _simulation;
	output_dir = std::move(_simulator_desc->output_dir);
	first_frame = _simulator_desc->first_frame;
	last_frame = _simulator_desc->last_frame;
	frames_per_second = _simulator_desc->frames_per_second;
	steps_per_frame = _simulator_desc->steps_per_frame;
	delete _simulator_desc;
}

void Simulator::Simulate()
{
	const double initial_time = omp_get_wtime();

	if (first_frame < 0) {
		Create_Output_Directory();
		Create_and_Write_Frame_Directory(0);

		std::string file_name = output_dir + "/scene_desc.txt";
		std::ofstream output(file_name);
		simulation->Write_Scene_Description(output);

		first_frame = 0;
	}
	else {
		Check_and_Load_Frame_Directory(first_frame);
	}
	
	// Initialize timing.
	const double begin_time = omp_get_wtime();
	std::cout << "#  Time: " << std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (begin_time - initial_time) << "s" << std::endl << std::endl;
	double last_time = begin_time;
	// Run the simulator.
	for (int frame = first_frame; frame < last_frame; frame++)
	{
		std::cout << "** Simulate Frame " << std::to_string(frame + 1) << "..." << std::endl;
		// Simulate.
		simulation->Reset_Time(Time_at_Frame(frame));
		Advance_Time_by_Steps(1.0 / frames_per_second);
		// Write files for frame.
		Create_and_Write_Frame_Directory(frame + 1);
		// Output timing.
		double current_time = omp_get_wtime();
		std::cout << "#  Time: "
			<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (current_time - last_time) << "s / "
			<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (current_time - initial_time) << "s" << std::endl;
		std::cout << "   Prediction:        "
			<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2)
			<< (current_time - begin_time) / (frame - first_frame + 1) * (last_frame - first_frame) + (begin_time - initial_time) << "s"
			<< std::endl << std::endl;
		last_time = current_time;
	}
}

void Simulator::Create_Output_Directory() const
{
	std::cout << "** Create output and write files for frame 0..." << std::endl;
	File::Create_Directory(output_dir);
}

void Simulator::Create_and_Write_Frame_Directory(const int frame) const
{
	// Create and write the frame directory.
	std::cout << "** Write output files for frame " << frame << "..." << std::endl;
	std::string frame_dir = output_dir + "/" + std::to_string(frame);
	File::Create_Directory(frame_dir);
	simulation->Write_Frame(frame_dir);
	simulation->Save_Frame(frame_dir);
	// Write the last frame.
	std::string file_name = output_dir + "/last_frame.txt";
	std::ofstream output(file_name);
	output << frame;
}

void Simulator::Check_and_Load_Frame_Directory(const int frame)
{
	std::cout << "** Load output files for frame " << frame << "..." << std::endl;
	std::string frame_dir = output_dir + "/" + std::to_string(frame);
	if (!File::Directory_Exists(frame_dir.c_str())) {
		std::cerr << "Error: [Simulator] No output directory for frame " << frame << "." << std::endl;
		std::exit(1);
	}
	simulation->Load_Frame(frame_dir);
}

void Simulator::Advance_Time_by_Steps(const double target_time)
{
	double time = 0.0;
	bool done = (time >= target_time);
	while (!done) {
		double begin_time = omp_get_wtime();

		double dt = Timestep();
		if (time + dt >= target_time) {
			dt = target_time - time;
			done = true;
		}
		else if (time + 2.0 * dt >= target_time) {
			dt = (target_time - time) * 0.5;
		}

		std::cout << "[" << std::setw(6) << static_cast<int>(target_time / dt + 0.5) << " SPF] ";
		simulation->Advance(dt);
		time += dt;
		double end_time = omp_get_wtime();
		std::cout << "    ... " << std::setw(7) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << end_time - begin_time << "s used" << std::endl;
	}
}
