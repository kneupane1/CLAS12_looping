
std::vector<float> *final_x_coord = new std::vector<float>();
std::vector<float> *final_y_coord = new std::vector<float>();

// you tringular cut ho, stafan ko rotation saga farak chha.

for (int i = 0; i < REC_Calorimeter_pindex->size(); i++) {
  if (REC_Calorimeter_pindex->size() == 0)
    continue;

  int j = REC_Calorimeter_pindex->at(i);
  if (REC_Particle_pid->at(j) == 11 && j == 0)

    if (REC_Calorimeter_sector->at(i) == 1 && REC_Calorimeter_x->at(i) < 350) {
      // printf(" Here I am\n");
      final_x_coord->push_back(REC_Calorimeter_x->at(i));
      final_y_coord->push_back(REC_Calorimeter_y->at(i));
    }

    else if (REC_Calorimeter_sector->at(i) == 2 &&
             REC_Calorimeter_x->at(i) > 0 &&
             REC_Calorimeter_y->at(i) > 0.577 * REC_Calorimeter_x->at(i) &&
             REC_Calorimeter_y->at(i) <
                 362 - (0.577 * REC_Calorimeter_x->at(i)))

    {
      final_x_coord->push_back(REC_Calorimeter_x->at(i));
      final_y_coord->push_back(REC_Calorimeter_y->at(i));

    } else if (REC_Calorimeter_sector->at(i) == 3 &&
               REC_Calorimeter_x->at(i) < -0.5 &&
               REC_Calorimeter_y->at(i) > -0.577 * REC_Calorimeter_x->at(i) &&
               REC_Calorimeter_y->at(i) <
                   362 + (0.577 * REC_Calorimeter_x->at(i)))

    {
      final_x_coord->push_back(REC_Calorimeter_x->at(i));
      final_y_coord->push_back(REC_Calorimeter_y->at(i));
    } else if (REC_Calorimeter_sector->at(i) == 4 &&
               REC_Calorimeter_x->at(i) > -362 &&
               REC_Calorimeter_y->at(i) <= -0.53 * (REC_Calorimeter_x->at(i))-44
}

{
  final_x_coord->push_back(REC_Calorimeter_x->at(i));
  final_y_coord->push_back(REC_Calorimeter_y->at(i));
}
else if (REC_Calorimeter_sector->at(i) == 5 && REC_Calorimeter_x->at(i) < 0 &&
         REC_Calorimeter_y->at(i) < 0.577 * REC_Calorimeter_x->at(i) &&
         REC_Calorimeter_y->at(i) > -362 - (0.577 * REC_Calorimeter_x->at(i)))

{
  final_x_coord->push_back(REC_Calorimeter_x->at(i));
  final_y_coord->push_back(REC_Calorimeter_y->at(i));
}
else if (REC_Calorimeter_sector->at(i) == 6 && REC_Calorimeter_x->at(i) > 0 &&
         REC_Calorimeter_y->at(i) < -0.577 * REC_Calorimeter_x->at(i) &&
         REC_Calorimeter_y->at(i) > -362 + (0.577 * REC_Calorimeter_x->at(i)))

{
  final_x_coord->push_back(REC_Calorimeter_x->at(i));
  final_y_coord->push_back(REC_Calorimeter_y->at(i));
}
}
for (int j = 0; j < final_x_coord->size(); j++) {
  pcal_Fiducial_cut_hist->Fill( // x_PCAL, y_PCAL);
      final_x_coord->at(j), final_y_coord->at(j));
  //  final_x_coord->clear();
  //  final_y_coord->clear();
}

/*
              if (REC_Calorimeter_sector->at(i) == 2) {
                if (REC_Calorimeter_x->at(i) > -362 &&
                        REC_Calorimeter_x->at(i) < 362 &&
                        REC_Calorimeter_y->at(i) > -362,
                    REC_Calorimeter_y->at(i) < 362 &&
                        TMath::Sqrt((REC_Calorimeter_x->at(i) *
                                     REC_Calorimeter_x->at(i)) +
                                    (REC_Calorimeter_y->at(i) *
                                     REC_Calorimeter_y->at(i))) <
   362) {  // yo mero circular cut thyo
                  */

// First particle is electron     (id[0] == ELECTRON);
// Number of good particles is greater than 0     (gpart > 0);
// First particle is negative     (q[0] == -1);
// First Particle is a good particle     (stat[0] > 0);
// First particle hit Electron Calorimeter    (ec[0] > 0);
// First Particle hit Scintillator     (sc[0] > 0);
// First Particle hit Drift Chamber     (dc[0] > 0);
// First Particle hit Cherenkov     (cc[0] > 0);

num_pip = 0;
good_e = false;
for (int j = 0; j < REC_Calorimeter_pindex->size(); j++) {
  if (REC_Calorimeter_pindex->size() == 0)
    continue;
  try {
    index = REC_Calorimeter_pindex->at(j);
    if (REC_Particle_pid->at(index) == 11) {
      e_mu_prime_3.SetXYZ(REC_Particle_px->at(index),
                          REC_Particle_py->at(index),
                          REC_Particle_pz->at(index));
      P = e_mu_prime_3.Mag();
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      sf = REC_Calorimeter_energy->at(j) / e_mu_prime.P();
      good_e = true;
    }
  } catch (std::exception &e) {
    total++;
  }
}
if (!good_e)
  continue;
good_e = false;
for (int j = 0; j < REC_Scintillator_time->size(); j++) {
  if (REC_Scintillator_time->size() == 0)
    continue;
  try {
    Delta_T *dt =
        new Delta_T(REC_Scintillator_time->at(0), REC_Scintillator_path->at(0));
    index = REC_Scintillator_pindex->at(j);
    sc_d = REC_Scintillator_detector->at(j);
    // I think 12 is FTOF
    if (sc_d == 12) {
      P_x = REC_Particle_px->at(index) * REC_Particle_px->at(index);
      P_y = REC_Particle_py->at(index) * REC_Particle_py->at(index);
      P_z = REC_Particle_pz->at(index) * REC_Particle_pz->at(index);
      P = TMath::Sqrt(P_x + P_y + P_z);

      dt->deltat(P, REC_Scintillator_time->at(j), REC_Scintillator_path->at(j));

      if (index == 0) {
        hist->Fill_MomVsBeta_vertex(REC_Particle_pid->at(index),
                                    REC_Particle_charge->at(index), P,
                                    beta->at(index));
        hist->Fill_deltat_vertex(REC_Particle_pid->at(index),
                                 REC_Particle_charge->at(index), P, dt);
      } else {
        hist->Fill_MomVsBeta(REC_Particle_pid->at(index),
                             REC_Particle_charge->at(index), P,
                             beta->at(index));
        hist->Fill_deltat(REC_Particle_pid->at(index),
                          REC_Particle_charge->at(index), P, dt);
      }
    }
    if (pid->at(REC_Scintillator_pindex->at(j)) == PIP &&
        abs(dt->Get_dt_Pi()) < 0.5)
      num_pip++;
    if (pid->at(REC_Scintillator_pindex->at(j)) == ELECTRON &&
        REC_Scintillator_detector->at(REC_Scintillator_pindex->at(j)) == 12)
      good_e = true;
    delete dt;
  } catch (std::exception &e) {
    total++;
  }
}
if (!good_e)
  continue;
// && sf >= 0.07 && sf <= 0.26
if (good_e && e_mu_prime.P() > 1.5) {
  Electromagnetic_Calorimeter->Fill(sf, e_mu_prime.P());
  W = physics::W_calc(e_mu, e_mu_prime);
  Q2 = physics::Q2_calc(e_mu, e_mu_prime);
  W_vs_q2->Fill(W, Q2);
  if (num_pip == 1 && pid->size() == 2)
    hist->Fill_WvsQ2_singlePi(W, Q2);
