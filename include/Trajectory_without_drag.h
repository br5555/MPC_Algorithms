#pragma once
#include <vector>
#include <math.h>

constexpr int num_time_segments = 8;
class Trajectory_without_drag
{
public:
	Trajectory_without_drag();
	Trajectory_without_drag(double jerk_max, double  jerk_min, double  acc_max, double  acc_min, double vel_max, double  vel_min, double  vel_start, double  vel_end, double  p_start, double  p_end
		, double  acc_start, double  acc_end);
	~Trajectory_without_drag();
private:
	bool calculate_trajectory_case1(double t[]);
	bool calculate_trajectory_case2(double t[]);
	bool calculate_trajectory_case3(double t[]);
	bool calculate_trajectory_case4(double t[]);
	bool calculate_trajectory_case5(double t[]);
	bool calculate_trajectory_case6(double t[]);
	bool calculate_trajectory_case7(double t[]);
	bool calculate_trajectory_case8(double t[]);
	std::vector<double> calculate_theta(double t[]);

	double jerk_max_, jerk_min_, acc_max_, acc_min_, vel_max_, vel_min_, vel_start_, vel_end_, p_start_, p_end_, acc_start_, acc_end_;

};

