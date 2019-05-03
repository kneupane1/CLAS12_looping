/*  if (REC_Particle_pid->at(j) == 11) {
    double px = REC_Particle_px->at(j) *
   REC_Particle_px->at(j);
    double py = REC_Particle_py->at(j) *
   REC_Particle_py->at(j);
    double pz = REC_Particle_pz->at(j) *
   REC_Particle_pz->at(j);
    P = TMath::Sqrt(px + py + pz);

    double energy = 0;
    for (int i = 0; i <
   REC_Calorimeter_pindex->size(); i++) {
      if (REC_Calorimeter_sector->at(i) ==
   1) {
        if (REC_Calorimeter_pindex->at(i)
   == 0)
   {
          energy = energy +
   REC_Calorimeter_energy->at(0);
        }
      }
    }

    SF = energy / P;
    if (SF > 0.05 && SF < 0.8) {
      SF_vs_P_hist->Fill(P, SF);
    }*/
// 2 decimal place mattry print out garna
cout.setf(ios::fixed);
cout.setf(ios::showpoint);
cout.precision(2);
