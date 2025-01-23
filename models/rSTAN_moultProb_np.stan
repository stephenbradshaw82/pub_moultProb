data {
  int<lower=1> nLobsters;                 //nLobsters;                              
  int<lower=1> nPeriods;                  //nPeriods (months);                              
  real obsGI[nLobsters];                  //Observation: Growth increment  
  int<lower=-1,upper=1> obsDam[nLobsters]; //Observation: Damage regeneration
  int<lower=-1,upper=1> obsPle[nLobsters]; //Observation: pleopod regeneration
  int<lower=-1,upper=1> obsMat[nLobsters]; //Observation: Maturity upgrade
  matrix<lower=0>[nLobsters,nPeriods] atLiberty; //time at liberty in each month
  int<lower=0> isFem;
}

parameters {
  real<lower=0, upper=1> p[nPeriods];
  real<lower=0, upper=100> growth_mu;
  real<lower=0.01, upper=10> growth_sigma;
  real<lower=0.5, upper=2> nogrowth_sigma;
  real<lower=0, upper=1> damageIndicates_moulted;
  real<lower=0, upper=1> damageIndicates_nonmoulted;
  real<lower=0, upper=1> pleopodIndicates_moulted;
  real<lower=0, upper=1> pleopodIndicates_nonmoulted;
  real<lower=0, upper=1> maturityIndicates_moulted;
  real<lower=0, upper=1> maturityIndicates_nonmoulted;
}

model {
  // //Weak Priors (Ziegler 2004, Macdiarmid 1989)
  // if(isFem==1){
  //     p[4] ~ beta(1.5, 1.5);
  //     p[5] ~ beta(1.5, 1.5);
  //     p[6] ~ beta(1.5, 1.5);
  //   } else {
  //     p[8] ~ beta(1.5, 1.5);
  //     p[9] ~ beta(1.5, 1.5);
  //     p[10] ~ beta(1.5, 1.5);
  // }//end isFem loop for p priors
  // 
  // //Informed Prior (Bradshaw 2024a)
  // nogrowth_sigma ~ normal(sqrt(2)*0.55, 0.05); 

  //Likelihood statement by lobster
  for (iLobster in 1:nLobsters){

    //Expression for moulting (time periods as months)
    real probNotMoulted = 1;
    real probmoulted;
    for (iPeriod in 1:nPeriods){
      probNotMoulted*=pow((1-p[iPeriod]),atLiberty[iLobster,iPeriod]);
    }
    probmoulted=1-probNotMoulted;

    // CONCEPT
    // (prob moulted) * p(growth incrment | moulted) *...
    // p(obs damage decrease | moulted) *...
    // p(obs pleopod regen | moulted) *...
    // p(obs maturity upgrade | moulted) ...
    // + 
    // (prob didn't moult) * p(growth increment|didn't moult) *...
    // prob(obs damage decrease | didn't moult) *...
    // p(obs pleopod regen | didn't moult) *...
    // p(obs maturity upgrade | didn't moult)
    
    //MOULTED & GROWTH INCREMENT EVIDENCE
    target += log(
      //Likelihood of moulting
      exp(
        bernoulli_lpmf(1 | probmoulted) +
        normal_lpdf(obsGI[iLobster] | growth_mu, growth_sigma)
      )+
      //Likelihood of NOT moulting
      exp(
        bernoulli_lpmf(0 | probmoulted) +
        normal_lpdf(obsGI[iLobster] | 0, nogrowth_sigma)
      )
    );

    //LIMB DAMAGE EVIDENCE
    if(obsDam[iLobster]>=0){
      target += log(
        //Likelihood of moulting --> LIMB REGEN
        exp(
          bernoulli_lpmf(1 | probmoulted) +
          bernoulli_lpmf(obsDam[iLobster] | damageIndicates_moulted)
        )+
        //Likelihood of NOT moulting --> LIMB REGEN
        exp(
          bernoulli_lpmf(0 | probmoulted) +
          bernoulli_lpmf(obsDam[iLobster] | damageIndicates_nonmoulted)
        )
      );
    }

    //PLEOPOD REGENERATION EVIDENCE
    if(obsPle[iLobster]>=0){
      target += log(
        //Likelihood of moulting --> PLEOPOD REGEN
        exp(
          bernoulli_lpmf(1 | probmoulted) +
          bernoulli_lpmf(obsPle[iLobster] | pleopodIndicates_moulted)
        )+
        //Likelihood of NOT moulting --> PLEOPOD REGEN
        exp(
          bernoulli_lpmf(0 | probmoulted) +
          bernoulli_lpmf(obsPle[iLobster] | pleopodIndicates_nonmoulted)
        )
      );
    }

    //MATURITY UPGRADE EVIDENCE
    if(obsMat[iLobster]>=0){
      target += log(
        //Likelihood of moulting --> MATURITY CHANGE 0->1
        exp(
          bernoulli_lpmf(1 | probmoulted) +
          bernoulli_lpmf(obsMat[iLobster] | maturityIndicates_moulted)
        )+
        //Likelihood of NOT moulting --> MATURITY CHANGE 0->1
        exp(
          bernoulli_lpmf(0 | probmoulted) +
          bernoulli_lpmf(obsMat[iLobster] | maturityIndicates_nonmoulted)
        )
      );
    }

  } //end for iLobster

} // end model block

