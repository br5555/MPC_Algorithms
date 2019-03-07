#include "LM_mpc_solver.h"



LM_MPC_Solver::LM_MPC_Solver(MatrixXd A, MatrixXd B, MatrixXd Bd, MatrixXd Q,
	MatrixXd Q_final, MatrixXd R, MatrixXd R_delta,
	MatrixXd disturbance, int num_params, int pred_horizon, int control_horizon, MatrixXd scale_MV, MatrixXd scale_OV) : num_params_(num_params), pred_horizon(pred_horizon), control_horizon(control_horizon), A_(A), B_(B), Bd_(Bd), Q_(Q), Q_final_(Q_final), R_(R), R_delta_(R_delta), disturbance_(disturbance)
{
	this->insecure_ = this->Bd_ * disturbance;
	x_states = MatrixXd::Zero(A_.rows(), pred_horizon + 1);
	deriv_wrt_u = MatrixXd::Zero(B_.cols(), control_horizon);
	u = MatrixXd::Zero(B_.cols(), 1);
	u_past = MatrixXd::Zero(B_.cols(), 1);
	u_current = MatrixXd::Zero(B_.cols(), 1);
	u_ss_ = MatrixXd::Zero(B_.cols(), 1);
	lambdas_x = MatrixXd::Zero(A_.rows(), pred_horizon);
	lambdas_u = MatrixXd::Zero(B_.cols(), pred_horizon);
	lambdas_u_ref = MatrixXd::Zero(B_.cols(), pred_horizon);
	u_horizon = MatrixXd::Zero(B_.cols(), pred_horizon);
	lambdas_u_ref = MatrixXd::Zero(B_.cols(), pred_horizon);
	x0_ = MatrixXd::Zero(x_states.rows(), 1);
	change_x = MatrixXd::Zero(x_states.rows(), 1);
	x_ss_ = MatrixXd::Zero(x_states.rows(), 1);
	Jacobian = Eigen::MatrixXd::Zero(1, 2 * control_horizon);
	mi_0 = 1.0;
	mi_dec = 0.1;
	mi_inc = 10.0;
	mi_max = 1e4;
	residuum = 0.0;
	max_iter = 30;
	alfa = 0.1;
	saturation_count = 0;
	min_residuum = 1e-5;
	scale_MV_inv = scale_MV.inverse();
	scale_OV_inv = scale_OV.inverse();
	A_pow_matrix = MatrixXd::Identity(A_.rows(), A_.cols());
	A_pow_B_cache = MatrixXd::Zero(A_.rows(), pred_horizon*B_.cols());
	A_pow_B_cache.block(0, 0, A_.rows(), B_.cols()) = MatrixXd::Identity(A_.rows(), A_.cols())* B_;
	/*//cout << "pow " << 0 << endl;
	//cout << A_pow_B_cache.block(0, (0), A_.rows(), B_.cols()) << endl;*/
	for (int i = 0; i < pred_horizon - 1; i++) {

		A_pow_B_cache.block(0, (i + 1)*B_.cols(), A_.rows(), B_.cols()) = (A_* A_pow_B_cache.block(0, (i)*B_.cols(), A_.rows(), B_.cols()));//.eval();
		/*//cout << "pow " << i + 1 << endl;
		//cout << A_pow_B_cache.block(0, (i + 1), A_.rows(), B_.cols()) << endl;*/
	}
	//cout << endl;
	////cout << A_pow_B_cache << endl;

}

//When you add the const keyword to a method the this pointer will essentially become a pointer to const object, and you cannot therefore change any member data. (This is not completely true, because you can mark a member as mutable and a const method can then change it. It's mostly used for internal counters and stuff.).

MatrixXd LM_MPC_Solver::Evaluate(Eigen::MatrixXd  x) {
	/*//cout << "R_delta " << endl << R_delta_ << endl;
	//cout << "R " << endl << R_ << endl;
	//cout << "Q " << endl << Q_ << endl;*/
	saturation_count = 0;
	mi_0 = 0.001;
	residuum_old = 1000.0;
	for (int iter = 0; iter < max_iter; iter++) {
		//cout << "Unutarnja iteracija " << iter <<endl;
		///TODO: Ove tri linije izbacit
		deriv_wrt_u.setZero();
		x_states.setZero();
		x_states.block(0, 0, x0_.rows(), x0_.cols()) = x0_;
		//cout << "X0 je "<<endl << x0_ << endl;
		//u_past.block(0,0,u_current.rows(), u_current.cols())= u_current;
		u_past = 0 * u_current; //  MatrixXd::Zero(B_.cols(), 1);


		/*//cout << "Derivacije" << endl;
		//cout << deriv_wrt_u << endl;
		//cout << "x_states" << endl <<  x_states << endl;
		//cout << "U_past " << endl << u_past << endl;*/
		//cout << "x_states" << endl << x_states << endl;

		for (int i = 0; i < pred_horizon; i++) {
			if (i < control_horizon) {
				u << x(0 * control_horizon + (i), 0), x(0 * control_horizon + (i), 0),
					x(1 * control_horizon + (i), 0), -x(1 * control_horizon + (i), 0);

			}
			//cout << "u " << endl << u << endl;
			/*//cout << "x " << endl << x << endl;*/

			x_states.block(0, i + 1, x0_.rows(), x0_.cols()) = (A_ * x_states.block(0, i, x0_.rows(), x0_.cols()) + B_ * u);//.eval();//+ insecure_;
			lambdas_x.block(0, i, x0_.rows(), x0_.cols()) = -1 * x_ss_ + x_states.block(0, i, x0_.rows(), x0_.cols());

			/*//cout << "x_states" << endl << x_states << endl;
			//cout << "lambdas_x" << endl << lambdas_x << endl;*/

			lambdas_u.block(0, i, u_past.rows(), u_past.cols()) = u - u_past;
			lambdas_u_ref.block(0, i, u.rows(), u.cols()) = u - u_ss_;

			/*out << "lambdas_u_ref" << endl << lambdas_u_ref << endl;*/

			//derivation of u
			if (i < control_horizon) {
				deriv_wrt_u.block(0, i, u.rows(), u.cols()) = (deriv_wrt_u.block(0, i, u.rows(), u.cols()) + (2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past));//.eval();
				//cout << "Derivacije prvi dio" << endl;
				//cout << deriv_wrt_u << endl;
				//cout << "Dodajem u prvi dio " << endl << (2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past) << endl;

				if (i > 0) {
					deriv_wrt_u.block(0, i - 1, u.rows(), u.cols()) = (deriv_wrt_u.block(0, i - 1, u.rows(), u.cols()) - 2 * R_delta_*u);//.eval();
					//cout << "Derivacije drugi dio" << endl;
					//cout << deriv_wrt_u << endl;
					//cout << "Dodajem u drugi dio " << endl << (- 2 * R_delta_*u) << endl;
				}
			}
			else {
				deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) = (deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) + (2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past));//.eval();
				//cout << "Derivacije treci dio" << endl;
				//cout << deriv_wrt_u << endl;
				//cout << "Dodajem u treci dio " << endl << (2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past) << endl;

				deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) = (deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) - 2 * R_delta_*u);//.eval();
				//cout << "Derivacije cetvrti dio" << endl;
				//cout << deriv_wrt_u << endl;
				//cout << "Dodajem u cetvrti dio " << endl << (deriv_wrt_u.block(0, 4, u.rows(), u.cols()) - 2 * R_delta_*u) << endl;

			}
			//derivation of x
			for (int j = 0; j <= i; j++) {


				if (j < control_horizon && i < control_horizon) {

					//bilo je promjenjen indeks u x_states_block
					//deriv_wrt_u.block(0, j, u.rows(), u.cols()) = (deriv_wrt_u.block(0, j, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j), A_.rows(), B_.cols())).transpose()).eval();
					deriv_wrt_u.block(0, j, u.rows(), u.cols()) = (deriv_wrt_u.block(0, j, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose());//.eval();
					/*//cout << i - j << endl;
					//cout << "Chache" << A_pow_B_cache.block(0, (i - j), A_.rows(), B_.cols()) << endl;*/
					//cout << "Derivacije peti dio" << endl;
					//cout << deriv_wrt_u << endl;
					//cout << "Dodajem u peti dio " << endl << ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose() << endl;
					//cout << "Peti dio prva zagrada" << endl << (2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()<< endl;
					//cout << "Peti dio prva zagrada" << endl << A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols()) << endl;
					//cout  << "x_states" << endl << x_states << endl;
					//cout << "samo dio" << endl;
					//cout << x_states.block(0, i + 1, x0_.rows(), x0_.cols()) << endl;
					//cout << "Samo provjera" << endl;
					//cout << x_ss_ << endl;
				}
				else {
					if (j >= control_horizon) {
						deriv_wrt_u.block(0, 4, u.rows(), u.cols()) = (deriv_wrt_u.block(0, 4, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose());//.eval();
						//cout << "Derivacije seti dio" << endl;
						//cout << deriv_wrt_u << endl;
						//cout << "Dodajem u sesti dio " << endl << ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose() << endl;
						
					}
					else {
						deriv_wrt_u.block(0, j, u.rows(), u.cols()) = (deriv_wrt_u.block(0, j, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose());//.eval();
						//cout << "Derivacije sedmi dio" << endl;
						//cout << deriv_wrt_u << endl;
						//cout << "Dodajem u sedmi dio " << endl << ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose() << endl;

					}

				}

			}

			//u_past.block(0,0,u_past.rows(), u_past.cols()) = u;
			u_past = u;
			//cout << "Derivacije" << endl;
			//cout << deriv_wrt_u << endl;
		}
		//cout << "Derivacije" << endl;
			//cout << deriv_wrt_u << endl;
		////cout << "jakobijan " << deriv_wrt_u << endl;
		lambdas_u_ref = scale_MV_inv * lambdas_u_ref;
		lambdas_u = scale_MV_inv * lambdas_u;
		lambdas_x = scale_OV_inv * lambdas_x;


		/*//cout << "lambdas_u_ref" << endl << lambdas_u_ref << endl;
		//cout << "lambdas_u" << endl << lambdas_u << endl;
		//cout << "lambdas_x" << endl << lambdas_x << endl;*/

		residuum_signal = (lambdas_u_ref.cwiseProduct(R_*lambdas_u_ref)).sum() + (lambdas_u.cwiseProduct(R_delta_*lambdas_u)).sum();

		residuum_state = (lambdas_x.cwiseProduct(Q_*lambdas_x)).sum();
		// +((-x_ss_ + x_states.block(0,horizon,x0_.rows(),x0_.cols())).transpose()*Q_final_*(-x_ss_ + x_states.block(0,horizon,x0_.rows(),x0_.cols()))).sum() );

		residuum = residuum_signal + residuum_state;//+ 1e1*(abs((lambdas_u_ref.block(0,0,1, pred_horizon) - lambdas_u_ref.block(1,0,1, pred_horizon)).sum())) +      1e1*(abs((lambdas_u_ref.block(2,0,1, pred_horizon) + lambdas_u_ref.block(3,0,1, pred_horizon)).sum()))             ;

		////cout << " residuum "<<residuum << "  "<<isnan(residuum) << endl;
		if (isnan(residuum)) residuum = 9e50;


		//ROS_INFO_STREAM("u_ss   : "<< u_ss_);
		//ROS_INFO_STREAM("x_ss_   : "<< x_ss_);
		//ROS_INFO_STREAM("residuum   : "<< residuum);

			//derivation df/du1
		Jacobian.block(0, 0, 1, control_horizon) = 2 * deriv_wrt_u.block(0, 0, 1, control_horizon);
		//derivation df/du2
		Jacobian.block(0, control_horizon, 1, control_horizon) = 2 * deriv_wrt_u.block(2, 0, 1, control_horizon);


		//cout << deriv_wrt_u << endl;
		//cout << endl;
		//cout << Jacobian << endl;
		//cout << endl;
		//cout << residuum << endl;
		//cout << endl;

		/*for (int iter_stab = 0; iter_stab < 2 * control_horizon; iter_stab++) {
			if (Jacobian(0, iter_stab) < 1e-10)
				Jacobian(0, iter_stab) = 0.0;
		}*/
		
		if (residuum < min_residuum) {
			return x;
		}



		



		//x = (x - 1 * (Jacobian.transpose() * Jacobian + mi_0 * Eigen::MatrixXd::Identity(2 * control_horizon, 2 * control_horizon)).inverse() * Jacobian.transpose() *residuum).eval();
		//Radi !!
		change_x = -0.1 * (Jacobian.transpose() * Jacobian + 10 * Eigen::MatrixXd::Identity(2 * control_horizon, 2 * control_horizon)).inverse() * Jacobian.transpose() *residuum;
		x = (x + change_x).eval();
		//cout << "X " << endl << x << endl;
		
		/*if (residuum <= residuum_old) {
			mi_0 = mi_0 * mi_dec;
		}
		else {
			mi_0 = mi_0 * mi_inc;

			if (mi_0 > mi_max) {
				return x;
			}
			else {
				iter--;
				continue;
			}

		}*/



		if (abs(residuum_old - residuum) < 1e-6) {
			saturation_count++;
			
		}
		else {
			saturation_count = 0;
		}

		if (saturation_count == 10) {
			return x;
		}

		if ((residuum -residuum_old  ) > 1e-4) {
			return x;
		}

		residuum_old = residuum;

	}


	return x;
}


LM_MPC_Solver::~LM_MPC_Solver()
{
}
