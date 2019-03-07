#include "Trajectory_without_drag.h"



Trajectory_without_drag::Trajectory_without_drag()
{
	double jerk_max_ = 5;
	double jerk_min_ = -5;
	double acc_max_ = 5.0;
	double acc_min_ = -5.0;
	double vel_max_ = 10;
	double vel_min_ = -10;
	double vel_start_ = 0;
	double vel_end_ = 0;
	double p_start_ = 0;
	double p_end_ = 10.3;
	double acc_start_ = 0;
	double acc_end_ = 0;
}

Trajectory_without_drag::Trajectory_without_drag(double jerk_max, double  jerk_min, double  acc_max, double  acc_min, double vel_max, double  vel_min, double  vel_start, double  vel_end, double  p_start, double  p_end
	, double  acc_start, double  acc_end) : jerk_max_(jerk_max), jerk_min_(jerk_min), acc_max_(acc_max), acc_min_(acc_min), vel_max_(vel_max), vel_min_(vel_min),
	vel_start_(vel_start), vel_end_(vel_end), p_start_(p_start), p_end_(p_end), acc_start_(acc_start), acc_end_(acc_end) {}

Trajectory_without_drag::~Trajectory_without_drag()
{
}

bool Trajectory_without_drag::calculate_trajectory_case1(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;

	t[0] = (acc_max_ - acc_start_) / jerk_max_;
	t[1] = (pow(acc_max_, 2.0)* jerk_max_ - pow(acc_max_, 2.0) * jerk_min_ + pow(acc_start_, 2.0) * jerk_min_ + 2.0 * jerk_min_* jerk_max_ * vel_max_ - 2.0 * jerk_min_ * jerk_max_ * vel_start_) / (2.0 * acc_max_ * jerk_min_ * jerk_max_);
	t[2] = -acc_max_ / jerk_min_;
	double A = acc_min_ * pow(acc_max_, 4.0) * pow(jerk_min_, 2.0) - pow(acc_min_, 4.0) * acc_max_ * pow(jerk_min_, 2.0) - acc_min_ * pow(acc_max_, 4.0) * pow(jerk_max_, 2.0) + pow(acc_min_, 4.0) * acc_max_*pow(jerk_max_, 2.0) - 3.0 * acc_min_ * pow(acc_start_, 4.0) * pow(jerk_min_, 2.0);
	double B = 3.0 * acc_max_ *pow(acc_end_, 4.0) * pow(jerk_min_, 2.0) + 8.0 * acc_min_ * acc_max_ * pow(acc_start_, 3.0) * pow(jerk_min_, 2.0) - 8.0 * acc_min_ * acc_max_ * pow(acc_end_, 3.0) * pow(jerk_min_, 2.0);

	double C = -6.0 * acc_min_ * pow(acc_max_, 2.0) * pow(acc_start_, 2.0) *pow(jerk_min_, 2.0) + 6.0 * pow(acc_min_, 2.0) * acc_max_ * pow(acc_end_, 2.0) *pow(jerk_min_, 2.0) + 12.0 * acc_min_ *pow(jerk_min_, 2.0) *pow(jerk_max_, 2.0) * pow(vel_max_, 2.0);

	double D = -12.0 * acc_max_ * pow(jerk_min_, 2.0) * pow(jerk_max_, 2.0) * pow(vel_max_, 2.0) - 12.0 * acc_min_ * pow(jerk_min_, 2.0) * pow(jerk_max_, 2.0) * pow(vel_start_, 2.0) + 12.0 * acc_max_ * pow(jerk_min_, 2.0) * pow(jerk_max_, 2.0) * pow(vel_end_, 2.0);

	double E = 24.0 * acc_min_ * acc_max_ * pow(jerk_min_, 2.0) * pow(jerk_max_, 2.0) * p_start_ - 24.0 * acc_min_ * acc_max_ * pow(jerk_min_, 2.0) * pow(jerk_max_, 2.0) * p_end_;

	double F = -12.0 * acc_min_ * pow(acc_max_, 2.0) * jerk_min_ * pow(jerk_max_, 2.0) * vel_max_ + 12.0 * pow(acc_min_, 2.0) * acc_max_ * jerk_min_ * pow(jerk_max_, 2.0) * vel_max_;

	double G = 12.0 * acc_min_ * pow(acc_max_, 2.0) *pow(jerk_min_, 2.0) * jerk_max_ * vel_start_ - 12.0 *pow(acc_min_, 2.0) * acc_max_ * pow(jerk_min_, 2.0) * jerk_max_ *vel_end_;

	double H = 12.0 * acc_min_ * pow(acc_start_, 2.0) * pow(jerk_min_, 2.0) * jerk_max_ * vel_start_ - 12.0 * acc_max_ *pow(acc_end_, 2.0) * pow(jerk_min_, 2.0) * jerk_max_ *vel_end_;

	double I = -24.0 * acc_min_ * acc_max_ * acc_start_ * pow(jerk_min_, 2.0) * jerk_max_ * vel_start_ + 24.0 * acc_min_ * acc_max_ * acc_end_ *pow(jerk_min_, 2.0) * jerk_max_ * vel_end_;

	double J = 24.0 * acc_min_ * acc_max_ * pow(jerk_min_, 2.0) * pow(jerk_max_, 2.0) * vel_max_;

	t[3] = -(A + B + C + D + E + F + G + H + I) / J;
	t[4] = acc_min_ / jerk_min_;
	t[5] = -(acc_min_, 2 * jerk_max_ - acc_min_, 2 * jerk_min_ + acc_end_, 2 * jerk_min_ + 2 * jerk_min_ * jerk_max_ * vel_max_ - 2 * jerk_min_ * jerk_max_ * vel_end_) / (2 * acc_min_ * jerk_min_ * jerk_max_);
	t[6] = -(acc_min_ - acc_end_) / (jerk_max_);

	for (int i = 0; i < num_time_segments; i++) {
		if (t[i] < 0)
			return calculate_trajectory_case2(t);
	}

	return true;
}

bool Trajectory_without_drag::calculate_trajectory_case2(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	
	/*
	double t1 = (acc_max_ - acc_start_) / jerk_max_;
	double a1 = acc_max_;
	double a5 = acc_min_;
	double t4 = 0;
	double a3 = 0;
	double a2 = acc_max_;
	double t3 = -acc_max_ / jerk_min_;
	double a4 = a3;
	double a7 = acc_end_;
	double a6 = acc_min_;
	double t7 = (a7 - a6) / jerk_max_;
	double a0 = acc_start_;
	double a5 = a6;
	double t5 = a5 / jerk_min_;
	double v1 = vel_start_ + t1 * acc_start_ + 0.5*pow(t1 , 2.0) * jerk_max_;
	double v6 = vel_end_ - t7 * a6 - 0.5*(pow((t7) , 2))*jerk_max_;
	double p6 = p_end_ - t7 * v6 - 0.5*(pow(t7 , 2.0))*a6 - (1.0 / 6.0)*(pow(t7 , 3.0))*jerk_max_;
	double p1 = p_start_ + t1 * vel_start_ + 0.5*(pow(t1 , 2.0))*acc_start_ + (1.0 / 6.0)*(pow(t1 , 3.0))*jerk_max_;
	double t6 = (2.0 * a1*v6 - 2.0 * a5*v6 + pow((4.0 * pow(a5 , 2.0) * pow(v1 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(v6 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * pow(t3 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * pow(t5 , 2.0) - 8.0 * a1*pow(a5 , 2.0) * p1 + 8.0 * pow(a1 , 2.0) * a5*p1 + 8.0 * a1*pow(a5 , 2.0) * p6 - 8.0 * pow(a1 , 2.0) * a5*p6 - 4.0 * a1*a5*pow(v1 , 2.0) - 4.0 * a1*a5*pow(v6 , 2.0) - 4.0 * a1*a2*pow(a5 , 2.0) * pow(t3 , 2.0) + 4.0 * a1*pow(a2 , 2.0) * a5*pow(t3 , 2.0) - 4.0 * pow(a1 , 2.0) * a2*a5*pow(t3 , 2.0) - (4.0 * a1*pow(a5 , 2.0) * jerk_min_*pow(t3 , 3.0)) / 3.0 - (8.0 * pow(a1 , 2.0) * a5*jerk_min_*pow(t3 , 3.0)) / 3.0 + a1 * a5*pow(jerk_min_ , 2.0) * pow(t3 , 4.0) - (4.0 * a1*pow(a5 , 2.0) * jerk_min_*pow(t5 , 3.0)) / 3.0 - (8.0 * pow(a1 , 2.0) * a5*jerk_min_*pow(t5 , 3.0)) / 3.0 + a1 * a5*pow(jerk_min_ , 2.0 )* pow(t5 , 4.0) + 8.0 *pow(a1 , 2.0) * pow(a5 , 2.0) * t3*t5 - 4.0 * a1*pow(a5 , 2.0) * jerk_min_*pow(t3 , 2.0) * t5 - 4.0 * pow(a1 , 2.0) * a5*jerk_min_*t3*pow(t5 , 2.0) + 2.0 * a1*a5*pow(jerk_min_ , 2.0) * pow(t3 , 2.0) * pow(t5 , 2.0) + 4.0 * a1*a2*a5*jerk_min_*pow(t3 , 3.0) - 8.0 * a1*a2*pow(a5 , 2.0) * t3*t5 + 4.0 * a1*a2*a5*jerk_min_*t3*pow(t5 , 2.0)) , (1.0 / 2.0)) - 2.0 * a1*a5*t3 + 2.0 * a2*a5*t3 - 2.0 * a1*a5*t5 + a5 * jerk_min_*pow(t3 , 2.0) + a5 * jerk_min_*pow(t5 , 2.0)) / (2.0 * a5*(a1 - a5));
	if (t6 < 0 )
		t6 = (2 * a1*v6 - 2 * a5*v6 - pow((4.0 * pow(a5 , 2.0) * pow(v1 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(v6 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * pow(t3 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * pow(t5 , 2.0) - 8.0 * a1*pow(a5 , 2.0) * p1 + 8.0 * pow(a1 , 2.0) * a5*p1 + 8.0 * a1*pow(a5 , 2.0) * p6 - 8.0 * pow(a1 , 2.0) * a5*p6 - 4.0 * a1*a5*pow(v1 , 2.0) - 4.0 * a1*a5*pow(v6 , 2.0) - 4.0 * a1*a2*pow(a5 , 2.0) * pow(t3 , 2.0) + 4.0 * a1*pow(a2 , 2.0) * a5*pow(t3 , 2.0) - 4.0 * pow(a1 , 2.0) * a2*a5*pow(t3 , 2.0) - (4.0 * a1*pow(a5 , 2.0) * jerk_min_*pow(t3 , 3.0)) / 3.0 - (8.0 * pow(a1 , 2.0) * a5*jerk_min_*pow(t3 , 3.0)) / 3.0 + a1 * a5*pow(jerk_min_ , 2.0) * pow(t3 , 4.0) - (4.0 * a1*pow(a5 , 2.0) * jerk_min_*pow(t5 , 3.0)) / 3.0 - (8.0 * pow(a1 , 2.0) * a5*jerk_min_*pow(t5 , 3.0)) / 3.0 + a1 * a5*pow(jerk_min_ , 2.0) * pow(t5 , 4.0) + 8.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * t3*t5 - 4.0 * a1*pow(a5 , 2.0) * jerk_min_*pow(t3 , 2.0) * t5 - 4.0 * pow(a1 , 2.0) * a5*jerk_min_*t3*pow(t5 , 2.0) + 2.0 * a1*a5*pow(jerk_min_ , 2.0) * pow(t3 , 2.0) * pow(t5 , 2.0) + 4.0 * a1*a2*a5*jerk_min_*pow(t3 , 3.0) - 8.0 * a1*a2*pow(a5 , 2.0) * t3*t5 + 4.0 * a1*a2*a5*jerk_min_*t3*pow(t5 , 2.0)) , (1.0 / 2.0)) - 2.0 * a1*a5*t3 + 2.0 * a2*a5*t3 - 2.0 * a1*a5*t5 + a5 * jerk_min_*pow(t3 , 2.0) + a5 * jerk_min_*pow(t5 , 2.0)) / (2.0 * a5*(a1 - a5));
	
		
	double 	v5 = v6 - t6 * a5;
	double v3 = v5 - 0.5*(pow(t5 , 2.0))*jerk_min_;
	double v4 = v3;
	double v2 = v3 - t3 * acc_max_ - 0.5*(pow(t3 , 2.0))*jerk_min_;

	double t2 = -(jerk_min_*pow(t3 , 2.0) + 2.0 * a2*t3 + jerk_min_ * pow(t5 , 2.0) + 2.0 * v1 - 2.0 * v6) / (2.0 * a1) - (2.0 * a1*v6 - 2.0 * a5*v6 + (4.0 * pow(a5 , 2.0) * pow(v1 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(v6 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * pow(t3 , 2.0) + 4.0 * pow(a1 , 2.0) * pow(a5 , 2.0) * pow(t5 , 2.0) - 8.0 * a1*pow(a5 , 2.0) * p1 + 8.0 * pow(a1 , 2.0) * a5*p1 + 8.0 * a1*pow(a5 , 2.0) * p6 - 8.0 * pow(a1 , 2.0) * a5*p6 - 4.0 * a1*a5*pow(v1 , 2.0) - 4.0 * a1*a5*pow(v6 , 2.0) - 4.0 * a1*a2*pow(a5 , 2.0) * pow(t3 , 2.0) + 4.0 * a1*pow(a2 , 2.0) * a5*pow(t3 , 2.0) - 4.0 * pow(a1 , 2.0) * a2*a5*t3 , 2 - (4 * a1*a5 , 2 * jerk_min_*t3 , 3) / 3 - (8 * a1 , 2 * a5*jerk_min_*t3 , 3) / 3 + a1 * a5*jerk_min_ , 2 * t3 , 4 - (4 * a1*a5 , 2 * jerk_min_*t5 , 3) / 3 - (8 * a1 , 2 * a5*jerk_min_*t5 , 3) / 3 + a1 * a5*jerk_min_ , 2 * t5 , 4 + 8 * a1 , 2 * a5 , 2 * t3*t5 - 4 * a1*a5 , 2 * jerk_min_*t3 , 2 * t5 - 4 * a1 , 2 * a5*jerk_min_*t3*t5 , 2 + 2 * a1*a5*jerk_min_ , 2 * t3 , 2 * t5 , 2 + 4 * a1*a2*a5*jerk_min_*t3 , 3 - 8 * a1*a2*a5 , 2 * t3*t5 + 4 * a1*a2*a5*jerk_min_*t3*t5 , 2) , (1 / 2) - 2 * a1*a5*t3 + 2 * a2*a5*t3 - 2 * a1*a5*t5 + a5 * jerk_min_*t3 , 2 + a5 * jerk_min_*t5 , 2) / (2 * a1*(a1 - a5));
	if (t2 < 0 )
		t2 = -(jerk_min_*t3 , 2 + 2 * a2*t3 + jerk_min_ * t5 , 2 + 2 * v1 - 2 * v6) / (2 * a1) - (2 * a1*v6 - 2 * a5*v6 - (4 * a5 , 2 * v1 , 2 + 4 * a1 , 2 * v6 , 2 + 4 * a1 , 2 * a5 , 2 * t3 , 2 + 4 * a1 , 2 * a5 , 2 * t5 , 2 - 8 * a1*a5 , 2 * p1 + 8 * a1 , 2 * a5*p1 + 8 * a1*a5 , 2 * p6 - 8 * a1 , 2 * a5*p6 - 4 * a1*a5*v1 , 2 - 4 * a1*a5*v6 , 2 - 4 * a1*a2*a5 , 2 * t3 , 2 + 4 * a1*a2 , 2 * a5*t3 , 2 - 4 * a1 , 2 * a2*a5*t3 , 2 - (4 * a1*a5 , 2 * jerk_min_*t3 , 3) / 3 - (8 * a1 , 2 * a5*jerk_min_*t3 , 3) / 3 + a1 * a5*jerk_min_ , 2 * t3 , 4 - (4 * a1*a5 , 2 * jerk_min_*t5 , 3) / 3 - (8 * a1 , 2 * a5*jerk_min_*t5 , 3) / 3 + a1 * a5*jerk_min_ , 2 * t5 , 4 + 8 * a1 , 2 * a5 , 2 * t3*t5 - 4 * a1*a5 , 2 * jerk_min_*t3 , 2 * t5 - 4 * a1 , 2 * a5*jerk_min_*t3*t5 , 2 + 2 * a1*a5*jerk_min_ , 2 * t3 , 2 * t5 , 2 + 4 * a1*a2*a5*jerk_min_*t3 , 3 - 8 * a1*a2*a5 , 2 * t3*t5 + 4 * a1*a2*a5*jerk_min_*t3*t5 , 2) , (1 / 2) - 2 * a1*a5*t3 + 2 * a2*a5*t3 - 2 * a1*a5*t5 + a5 * jerk_min_*t3 , 2 + a5 * jerk_min_*t5 , 2) / (2 * a1*(a1 - a5));
	
		
	//t2 = (v2 - v1) / a1;

	*/
	return true;

}

bool Trajectory_without_drag::calculate_trajectory_case3(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	return true;

}
bool Trajectory_without_drag::calculate_trajectory_case4(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	return true;

}
bool Trajectory_without_drag::calculate_trajectory_case5(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	return true;

}
bool Trajectory_without_drag::calculate_trajectory_case6(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	return true;
}
bool Trajectory_without_drag::calculate_trajectory_case7(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	return true;
}
bool Trajectory_without_drag::calculate_trajectory_case8(double t[]) {
	double times_[num_time_segments] = {};
	t = times_;
	return true;
}