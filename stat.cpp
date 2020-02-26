  #include "stat.hpp"
  
  Stat::Stat(void){
      summ = 0.0f;
      sq_summ = 0.0f;
      var = 0.0f;
      std = 0.0f;
      summ_0 = 0.0f;
      summ_xi = 0.0f;
      sqsumm_xi = 0.0f;
      var_xi = 0.0f;
      crossumm = 0.0f;
      covar = 0.0f;
      counter = 0;
      counterN  = 0;
      summN = 0;
  }

  float Stat::Update_FKAK(Eigen::VectorXf X, float Y, float Z, float xi, int N){
      summ +=
  }

  float Stat::Fast_update(Eigen::VectorXf X, float Y, float Z, float xi){



  }