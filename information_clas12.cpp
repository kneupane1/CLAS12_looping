
#ifndef MAIN_H_GUARD // kuny kura 2 choti define nabhawos bhanera.
#define MAIN_H_GUARD
#include "TChain.h"
#include "TF1.h"
#include "TH2.h"
#include "deltat.hpp"
#include <TFile.h>
#include <TLorentzVector.h>
#include <fstream>
#include <vector> // vector is dynamic array, unlike array this can be resize. we can add elements: start from 1 and make 2,3 or 10 or ...;

std::vector<int> *REC_Particle_pid;
std::vector<float> *REC_Particle_px;
std::vector<float> *REC_Particle_py;
std::vector<float> *REC_Particle_pz;
std::vector<float> *REC_Particle_vx;
std::vector<float> *REC_Particle_vy;
std::vector<float> *REC_Particle_vz;
std::vector<int> *REC_Particle_charge;
std::vector<float> *REC_Particle_beta;
std::vector<float> *REC_Particle_chi2pid;
std::vector<int> *REC_Particle_status;

std::vector<int> *REC_Calorimeter_index;
std::vector<int> *REC_Calorimeter_pindex;
std::vector<int> *REC_Calorimeter_detector;
std::vector<int> *REC_Calorimeter_sector;
std::vector<int> *REC_Calorimeter_layer;
std::vector<float> *REC_Calorimeter_energy;
std::vector<float> *REC_Calorimeter_time;
std::vector<float> *REC_Calorimeter_path;
std::vector<float> *REC_Calorimeter_chi2;
std::vector<float> *REC_Calorimeter_x;
std::vector<float> *REC_Calorimeter_y;
std::vector<float> *REC_Calorimeter_z;
std::vector<float> *REC_Calorimeter_hx;
std::vector<float> *REC_Calorimeter_hy;
std::vector<float> *REC_Calorimeter_hz;
std::vector<float> *REC_Calorimeter_lu;
std::vector<float> *REC_Calorimeter_lv;
std::vector<float> *REC_Calorimeter_lw;
std::vector<float> *REC_Calorimeter_du;
std::vector<float> *REC_Calorimeter_dv;
std::vector<float> *REC_Calorimeter_dw;
std::vector<float> *REC_Calorimeter_m2u;
std::vector<float> *REC_Calorimeter_m2v;
std::vector<float> *REC_Calorimeter_m2w;
// std::vector<float> *REC_Calorimeter_m3u;
// std::vector<float> *REC_Calorimeter_m3v;
// std::vector<float> *REC_Calorimeter_m3w;
std::vector<int> *REC_Calorimeter_status;

std::vector<int> *REC_Scintillator_status;
std::vector<int> *REC_Scintillator_pindex;
std::vector<int> *REC_Scintillator_index;
std::vector<float> *REC_Scintillator_time;
std::vector<float> *REC_Scintillator_path;

std::vector<int> *REC_Cherenkov_index;
std::vector<int> *REC_Cherenkov_pindex;
std::vector<int> *REC_Cherenkov_status;
std::vector<float> *REC_Cherenkov_time;
std::vector<float> *REC_Cherenkov_path;

std::vector<int> *REC_ForwardTagger_index;
std::vector<int> *REC_ForwardTagger_pindex;
std::vector<int> *REC_ForwardTagger_status;
std::vector<float> *REC_ForwardTagger_time;
std::vector<float> *REC_ForwardTagger_path;

void test(char *fin, char *fout) {
  TFile *out = new TFile(fout, "RECREATE");
  double P;
  double SF;
  bool electron_cuts;

  // Load chain from branch h10
  TChain chain("clas12");
  chain.Add(fin);

  chain.SetBranchAddress("REC_Particle_pid", &REC_Particle_pid);
  chain.SetBranchAddress("REC_Particle_px", &REC_Particle_px);
  chain.SetBranchAddress("REC_Particle_py", &REC_Particle_py);
  chain.SetBranchAddress("REC_Particle_pz", &REC_Particle_pz);
  chain.SetBranchAddress("REC_Particle_vx", &REC_Particle_vx);
  chain.SetBranchAddress("REC_Particle_vy", &REC_Particle_vy);
  chain.SetBranchAddress("REC_Particle_vz", &REC_Particle_vz);
  chain.SetBranchAddress("REC_Particle_charge", &REC_Particle_charge);
  chain.SetBranchAddress("REC_Particle_beta", &REC_Particle_beta);
  chain.SetBranchAddress("REC_Particle_chi2pid", &REC_Particle_chi2pid);
  chain.SetBranchAddress("REC_Particle_status", &REC_Particle_status);

  chain.SetBranchAddress("REC_Calorimeter_energy", &REC_Calorimeter_energy);
  chain.SetBranchAddress("REC_Calorimeter_index", &REC_Calorimeter_index);
  chain.SetBranchAddress("REC_Calorimeter_status", &REC_Calorimeter_status);
  chain.SetBranchAddress("REC_Calorimeter_pindex", &REC_Calorimeter_pindex);
  chain.SetBranchAddress("REC_Calorimeter_detector", &REC_Calorimeter_detector);
  chain.SetBranchAddress("REC_Calorimeter_sector", &REC_Calorimeter_sector);
  chain.SetBranchAddress("REC_Calorimeter_layer", &REC_Calorimeter_layer);
  chain.SetBranchAddress("REC_Calorimeter_time", &REC_Calorimeter_time);
  chain.SetBranchAddress("REC_Calorimeter_path", &REC_Calorimeter_path);
  chain.SetBranchAddress("REC_Calorimeter_chi2", &REC_Calorimeter_chi2);
  chain.SetBranchAddress("REC_Calorimeter_x", &REC_Calorimeter_x);
  chain.SetBranchAddress("REC_Calorimeter_y", &REC_Calorimeter_y);
  chain.SetBranchAddress("REC_Calorimeter_z", &REC_Calorimeter_z);
  chain.SetBranchAddress("REC_Calorimeter_hx", &REC_Calorimeter_hx);
  chain.SetBranchAddress("REC_Calorimeter_hy", &REC_Calorimeter_hy);
  chain.SetBranchAddress("REC_Calorimeter_hz", &REC_Calorimeter_hz);
  chain.SetBranchAddress("REC_Calorimeter_lu", &REC_Calorimeter_lu);
  chain.SetBranchAddress("REC_Calorimeter_lv", &REC_Calorimeter_lv);
  chain.SetBranchAddress("REC_Calorimeter_lw", &REC_Calorimeter_lw);
  chain.SetBranchAddress("REC_Calorimeter_du", &REC_Calorimeter_du);
  chain.SetBranchAddress("REC_Calorimeter_dv", &REC_Calorimeter_dv);
  chain.SetBranchAddress("REC_Calorimeter_dw", &REC_Calorimeter_dw);
  chain.SetBranchAddress("REC_Calorimeter_m2u", &REC_Calorimeter_m2u);
  chain.SetBranchAddress("REC_Calorimeter_m2v", &REC_Calorimeter_m2v);
  chain.SetBranchAddress("REC_Calorimeter_m2w", &REC_Calorimeter_m2w);
  // chain.SetBranchAddress("REC_Calorimeter_m3u", &REC_Calorimeter_m3u);
  // chain.SetBranchAddress("REC_Calorimeter_m3v", &REC_Calorimeter_m3v);
  // chain.SetBranchAddress("REC_Calorimeter_m3w", &REC_Calorimeter_m3w);

  chain.SetBranchAddress("REC_Scintillator_pindex", &REC_Scintillator_pindex);
  chain.SetBranchAddress("REC_Scintillator_status", &REC_Scintillator_status);
  chain.SetBranchAddress("REC_Scintillator_path", &REC_Scintillator_path);
  chain.SetBranchAddress("REC_Scintillator_time", &REC_Scintillator_time);
  chain.SetBranchAddress("REC_Scintillator_index", &REC_Scintillator_index);
  chain.SetBranchAddress("REC_Cherenkov_pindex", &REC_Cherenkov_pindex);
  chain.SetBranchAddress("REC_Cherenkov_index", &REC_Cherenkov_index);
  chain.SetBranchAddress("REC_Cherenkov_status", &REC_Cherenkov_status);
  chain.SetBranchAddress("REC_Cherenkov_path", &REC_Cherenkov_path);
  chain.SetBranchAddress("REC_Cherenkov_time", &REC_Cherenkov_time);
  chain.SetBranchAddress("REC_ForwardTagger_pindex", &REC_ForwardTagger_pindex);
  chain.SetBranchAddress("REC_ForwardTagger_index", &REC_ForwardTagger_index);
  chain.SetBranchAddress("REC_ForwardTagger_status", &REC_ForwardTagger_status);
  chain.SetBranchAddress("REC_ForwardTagger_time", &REC_ForwardTagger_time);
  chain.SetBranchAddress("REC_ForwardTagger_path", &REC_ForwardTagger_path);

  int num_of_events = (int)chain.GetEntries();

  int total = 0;
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    if (REC_Particle_pid->size() == 0 || REC_Calorimeter_pindex->size() == 0 ||
        REC_Scintillator_time->size() == 0 || REC_Cherenkov_pindex->size() == 0)
      continue;

    double per = ((double)current_event / (double)num_of_events);

    std::cerr << "\t\t" << current_event << "\t\t"
              << std::floor(
                     (100 * (double)current_event / (double)num_of_events))
              << "%\r\r" << std::flush;

    for (int i = 0; i < 17; i++) {
      std::cout << REC_Particle_pid->at(i) << "      " << REC_Particle_px->at(i)
                << "      " << REC_Particle_py->at(i) << "      "
                << REC_Particle_pz->at(i) << "      " << REC_Particle_vx->at(i)
                << "      " << REC_Particle_vy->at(i) << "      "
                << REC_Particle_vz->at(i) << "      "
                << REC_Particle_charge->at(i) << "      "
                << REC_Particle_beta->at(i) << "      "
                << REC_Particle_chi2pid->at(i) << "      "
                << REC_Particle_status->at(i) << "      " <<
    }

    for (int i = 0; i < 17; i++) {
      std::cout
          << REC_Calorimeter_index->at(i) << "   "
          << REC_Calorimeter_pindex->at(i) << "   "
          << REC_Calorimeter_detector->at(i) << "   "
          << REC_Calorimeter_sector->at(i) << "   "
          << REC_Calorimeter_layer->at(i) << "   "
          << REC_Calorimeter_energy->at(i) << "   "
          << REC_Calorimeter_time->at(i) << "   " << REC_Calorimeter_path->at(i)
          << "   " << REC_Calorimeter_chi2->at(i) << "   "
          << REC_Calorimeter_x->at(i) << "    " << REC_Calorimeter_y->at(i)
          << "     " << REC_Calorimeter_z->at(i) << "    "
          << REC_Calorimeter_hx->at(i) << "   " << REC_Calorimeter_hy->at(i)
          << "   " << REC_Calorimeter_hz->at(i) << "   "
          << REC_Calorimeter_lu->at(i) << "   " << REC_Calorimeter_lv->at(i)
          << "   " << REC_Calorimeter_lw->at(i) << "   "
          << REC_Calorimeter_du->at(i) << "   " << REC_Calorimeter_dv->at(i)
          << "   " << REC_Calorimeter_dw->at(i) << "   "
          << REC_Calorimeter_m2u->at(i) << "   " << REC_Calorimeter_m2v->at(i)
          << "   " << REC_Calorimeter_m2w->at(i) << "   "
          //<< REC_Calorimeter_m3u->at(i) << "   " << REC_Calorimeter_m3v->at(i)
          //  << "   " << REC_Calorimeter_m3w->at(i) << "   "
          << REC_Calorimeter_status->at(i) <<

          '\n';
    }

    /* << "  " << REC_Particle_py->at(i) << "  "
    << REC_Particle_pz->at(i) << "  " << REC_Particle_vx->at(i)
    << "  " << REC_Particle_charge->at(i) << "  "
    << REC_Particle_beta->at(i) << "  "
    << REC_Particle_status->at(i) << "  "
    /*<< REC_Calorimeter_energy->at(i) << "  "
    << REC_Calorimeter_pindex->at(i) << "  "
    << REC_Scintillator_pindex->at(i) << "  "
    << REC_Scintillator_path->at(i) << "  "
    << REC_Cherenkov_pindex->at(i) << "  "
    << REC_Cherenkov_path->at(i) << "  "

        <<

            '\n';
      } */

    chain.Reset();
    std::cerr << "\n" << total << std::endl;
  }
#endif
