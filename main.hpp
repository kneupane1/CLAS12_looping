
#ifndef MAIN_H_GUARD // kuny kura 2 choti define nabhawos bhanera.
#define MAIN_H_GUARD
#include "TChain.h"
#include "TF1.h"
#include "TH2.h"
#include "colors.hpp"
#include "deltat.hpp"
#include <TFile.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <cmath>
#include <cstdlib>
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
std::vector<int> *REC_Scintillator_detector;

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
std::vector<float> *REC_ForwardTagger_x;
std::vector<float> *REC_ForwardTagger_y;

std::vector<int> *REC_Track_pindex;
std::vector<int> *REC_Track_detector;
std::vector<int> *REC_Track_sector;

std::vector<int> *REC_Traj_pindex;
std::vector<int> *REC_Traj_detId;
std::vector<float> *REC_Traj_x;
std::vector<float> *REC_Traj_y;
std::vector<float> *REC_Traj_z;

static const double Pival = TMath::Pi();

TH1D *momentum = new TH1D("Momentum", "Momentum; P (GeV); Entries", 500, 0,
                          2.40); //(const char *name, const char *title, Int_t
// nbinsx, Double_t xlow, Double_t xup)

TH1D *W_hist = new TH1D("W", "W ; W (GeV); Entries", 500, 0.8, 1.8);

TH1D *Q2_hist = new TH1D("Q2", "Q2; Q2 (GeV2); Entries", 500, 0.0, 1.0);
TH2D *W_vs_q2 =
    new TH2D("W_vs_Q2", "W_vs_Q2 ; W (GeV); Q2 (GeV2)", 500, 0.8, 1.9, 500, 0,
             1); //(const char *name, const char *title, Int_t
                 // nbinsx, Double_t xlow, Double_t xup, Int_t
// nbinsy, Double_t ylow, Double_t yup)
// TH2D *mom_vs_beta_0 = new TH2D("mom_vs_beta_0", "mom_vs_beta_0; P (GeV);
// beta",
// 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_except_e = new TH2D(
    "mom_vs_beta", "mom_vs_beta;  P (GeV);beta ", 500, 0, 2.5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_pos =
    new TH2D("mom_vs_beta_pos", "mom_vs_beta_pos; P (GeV); beta+ ", 500, 0, 2.5,
             500, 0.0, 1.2);
TH2D *mom_vs_beta_neg =
    new TH2D("mom_vs_beta_neg", "mom_vs_beta_neg; P (GeV); beta-", 500, 0, 2.5,
             500, 0.0, 1.2);
TH2D *mom_vs_beta_proton =
    new TH2D("mom_vs_beta_proton", "mom_vs_beta_proton; P (GeV); beta ", 500, 0,
             2.5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_anti_proton = new TH2D(
    "mom_vs_beta_anti_proton", "mom_vs_beta_anti_proton; P (GeV); beta ", 500,
    0, 2.5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_pion =
    new TH2D("mom_vs_beta_pion", "mom_vs_beta_pion; P (GeV); beta", 500, 0, 2.5,
             500, 0.0, 1.2);
TH2D *mom_vs_beta_pion_mi =
    new TH2D("mom_vs_beta_pion_mi", "mom_vs_beta_pion_mi; P (GeV); beta", 500,
             0, 2.5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_electron =
    new TH2D("mom_vs_beta_electron", "mom_vs_beta_electron; P (GeV); beta", 500,
             0, 2.5, 500, 0.0, 1.2);

TH1D *PCAL_hist = new TH1D(" PCAL", "PCAL ; PCAL (GeV); Entries", 500, 0, 0.6);
TH1D *Z_Vertex_position_hist = new TH1D(
    " Z-Vertex position", "Z-Vertex position ; Vz(/cm); Entries", 500, -40, 50);

TH1D *ECAL_hist = new TH1D(" ECAL", "ECAL ; ECAL (GeV); Entries", 500, 0, 0.6);
TH2D *SF_vs_P_hist = new TH2D("SF_vs_P", "SF_vs_P; P (GeV);	SF ", 500, 0.8,
                              2.3, 500, 0.0, 0.55);
TH2D *pcal_SF_vs_P_hist =
    new TH2D("pcal_SF_vs_P", "pcal_SF_vs_P; P (GeV);	pcal_SF ", 500, 0.6,
             2.3, 500, 0.0, 0.55);
TH2D *ecal_SF_vs_P_hist =
    new TH2D("ecal_SF_vs_P", "ecal_SF_vs_P; P (GeV);	ecal_SF ", 500, 0.6,
             2.3, 500, 0.0, 0.55);
TH2D *ECAL_vs_PCAL_hist =
    new TH2D("ECAL_vs_PCAL", "ECAL_vs_PCAL; PCAL (GeV);	ECAL (GeV)", 500, 0,
             0.6, 500, 0, 0.6);
TH2D *ECout_vs_ECin_hist =
    new TH2D("ECout_vs_ECin", "ECout_vs_ECin; ECin (GeV);	ECout (Gev) ",
             500, 0, 0.5, 500, 0, 0.5);
TH2D *pcal_Fiducial_cut_hist =
    new TH2D("pcal_Fiducial_cut", "pcal_Fiducial_cut; x/cm;	y/cm ", 1500,
             -400, 400, 1500, -400, 400);

TH2D *dc_Fiducial_cut_hist =
    new TH2D("dc_Fiducial_cut", "dc_Fiducial_cut; x/cm;	y/cm ", 1500, -400, 400,
             1500, -400, 400);

TH2D *delta_t_electron_vertex_withoutID_vs_P_hist =
    new TH2D("delta_t_electron_vertex__withoutID_vs_P",
             "delta_t_electron_vertex_withoutID_vs_P; P (GeV); delta_t (ns)",
             500, 0.0, 2.2, 500, -10, 10);

TH2D *delta_t_electron_vertex_withID_vs_P_hist =
    new TH2D("delta_t_electron_vertex_withID_vs_P",
             "delta_t_electron_vertex_withID_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -10, 10);

TH2D *delta_t_electron_vertex_antiID_vs_P_hist =
    new TH2D("delta_t_electron_vertex_antiID_vs_P",
             "delta_t_electron_vertex_antiID_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -10, 10);

TH2D *delta_t_electron_vs_P_hist = new TH2D(
    "delta_t_electron_vs_P", "delta_t_electron_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -20, 20);

TH2D *delta_t_anti_electron_vs_P_hist =
    new TH2D("delta_t_electron_anti_vs_P",
             "delta_t_anti_electron_vs_P; P (GeV); delta_t (ns)", 500, 0.0, 2.2,
             500, -20, 20);
TH2D *delta_t_proton_vs_P_hist = new TH2D(
    "delta_t_proton_vs_P", "delta_t_proton_vs_P; P (GeV); delta_t (ns)", 500,
    0.0, 2.2, 500, -20, 20);

TH1D *proton_P_hist = new TH1D(
    " Momentum_p", "Momentum_p ; Momentum_p (GeV); Entries", 500, 0, 2);

TH2D *delta_t_anti_proton_vs_P_hist =
    new TH2D("delta_t_anti_proton_vs_P",
             "delta_t_anti_proton_vs_P; P (GeV); delta_t (ns)", 500, 0.0, 2.2,
             500, -30, 30);

TH2D *delta_t_anti_pip_vs_P_hist = new TH2D(
    "delta_t_anti_pip_vs_P", "delta_t_anti_pip_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -30, 30);

TH2D *delta_t_pip_vs_P_hist =
    new TH2D("delta_t_pip_vs_P", "delta_t_pip_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -30, 30);
TH2D *delta_t_pi0_vs_P_hist =
    new TH2D("delta_t_pi0_vs_P", "delta_t_pi0_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -20, 20);

TH2D *delta_t_pim_vs_P_hist =
    new TH2D("delta_t_pim_vs_P", "delta_t_pim_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -20, 20);
TH2D *delta_t_kap_vs_P_hist =
    new TH2D("delta_t_kap_vs_P", "delta_t_kap_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -30, 30);
TH2D *delta_t_kam_vs_P_hist =
    new TH2D("delta_t_kam_vs_P", "delta_t_kam_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -20, 20);

TH2D *delta_t_g_vs_P_hist =
    new TH2D("delta_t_g_vs_P", "delta_t_g_vs_P; P (GeV); delta_t (ns)", 500,
             0.0, 2.2, 500, -10, 10);
TH2D *delta_t_omega_vs_P_hist =
    new TH2D("delta_t_omega_vs_P", "delta_t_omega_vs_P; P (GeV); delta_t (ns)",
             500, 0.0, 2.2, 500, -10, 10);

TH2D *delta_t_electron_vs_P_Withid_hist =
    new TH2D("delta_t_electron_Withid_vs_P",
             "delta_t_electron_vs_P_Withid; P (GeV); delta_t (ns)", 500, 0.0,
             2.2, 500, -20, 20);

TH2D *delta_t_proton_vs_P_Withid_hist =
    new TH2D("delta_t_proton_Withid_vs_P",
             "delta_t_proton_Withid_vs_P; P (GeV); delta_t (ns)", 500, 0.0, 2.2,
             500, -20, 20);

TH2D *delta_t_pip_vs_P_Withid_hist = new TH2D(
    "delta_t_pip_Withid_vs_P", "delta_t_pip_withid_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -30, 30);
TH2D *delta_t_pi0_vs_P_Withid_hist = new TH2D(
    "delta_t_pi0_Withid_vs_P", "delta_t_pi0_Withid_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -10, 10);

TH2D *delta_t_pim_vs_P_Withid_hist = new TH2D(
    "delta_t_pim_Withid_vs_P", "delta_t_pim_Withid_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -10, 10);
TH2D *delta_t_kap_vs_P_Withid_hist = new TH2D(
    "delta_t_kap_Withid_vs_P", "delta_t_kap_Withid_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -30, 30);
TH2D *delta_t_kam_vs_P_Withid_hist = new TH2D(
    "delta_t_kam_Withid_vs_P", "delta_t_kam_Withid_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -20, 20);

TH2D *delta_t_g_vs_P_Withid_hist = new TH2D(
    "delta_t_g_Withid_vs_P", "delta_t_g_Withid_vs_P; P (GeV); delta_t (ns)",
    500, 0.0, 2.2, 500, -10, 10);
// Calcuating Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}

//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  p_mu.SetVectM(p_mu_3, MASS_P);
  return (p_mu + q_mu).Mag();
}
void test(char *fin, char *fout) {
  TFile *out = new TFile(fout, "RECREATE");
  double P;
  double SF;
  double pcal_SF;
  double ecal_SF;
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

  chain.SetBranchAddress("REC_Traj_pindex", &REC_Traj_pindex);
  chain.SetBranchAddress("REC_Traj_detId", &REC_Traj_detId);
  chain.SetBranchAddress("REC_Traj_x", &REC_Traj_x);
  chain.SetBranchAddress("REC_Traj_y", &REC_Traj_y);
  chain.SetBranchAddress("REC_Traj_z", &REC_Traj_z);

  chain.SetBranchAddress("REC_Track_pindex", &REC_Track_pindex);
  chain.SetBranchAddress("REC_Track_detector", &REC_Track_detector);
  chain.SetBranchAddress("REC_Track_sector", &REC_Track_sector);

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
  //  chain.SetBranchAddress("REC_Calorimeter_m2w", &REC_Calorimeter_m2w);
  // chain.SetBranchAddress("REC_Calorimeter_m3u", &REC_Calorimeter_m3u);
  // chain.SetBranchAddress("REC_Calorimeter_m3v", &REC_Calorimeter_m3v);
  // chain.SetBranchAddress("REC_Calorimeter_m3w", &REC_Calorimeter_m3w);

  chain.SetBranchAddress("REC_Scintillator_pindex", &REC_Scintillator_pindex);
  chain.SetBranchAddress("REC_Scintillator_status", &REC_Scintillator_status);
  chain.SetBranchAddress("REC_Scintillator_path", &REC_Scintillator_path);
  chain.SetBranchAddress("REC_Scintillator_time", &REC_Scintillator_time);
  chain.SetBranchAddress("REC_Scintillator_index", &REC_Scintillator_index);
  chain.SetBranchAddress("REC_Scintillator_detector",
                         &REC_Scintillator_detector);
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
  chain.SetBranchAddress("REC_ForwardTagger_x", &REC_ForwardTagger_x);
  chain.SetBranchAddress("REC_ForwardTagger_y", &REC_ForwardTagger_y);

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

    for (int j = 0; j < REC_Particle_pid->size(); j++) {
      if (REC_Particle_pid->at(j) == 11) {
        //      std::cout << "my det dc " << REC_Traj_x->at(j);
        if (REC_Particle_charge->at(0) == -1) {
          if (REC_Particle_vz->at(j) != 0) {
            Z_Vertex_position_hist->Fill(REC_Particle_vz->at(j));
          }
        }
      }
    }
    // if (REC_Particle_pid->at(j) == 11 && j == 0)

    //  if (REC_Scintillator_pindex->at(index) == 11) {

    //          for (int k = 0; k < REC_Scintillator_pindex->size();
    // k++)
    //          {
    //          if (REC_Scintillator_pindex->at(k) == 11) {

    /*  for (int l = 0; l < REC_Cherenkov_pindex->size(); l++) {
        if (REC_Cherenkov_pindex->at(l) == ) {

          for (int m = 0; m < REC_ForwardTagger_pindex->size(); m++) {
            if (REC_ForwardTagger_pindex->at(m) == ) {
              /*if (REC_ForwardTagger_index->size() == 0 ||
                  REC_Calorimeter_index->size() == 0 ||
                  REC_Scintillator_index->size() == 0 ||
                  REC_Cherenkov_index->size() == 0)
                continue;
              /*  if (REC_Scintillator_status->at(0) == 11) {

                  if (REC_Cherenkov_status->at(0) == 11) {

                    if (REC_ForwardTagger_status->at(0) ==
              11) {*/

    for (int i = 1; i < REC_Particle_pid->size(); i++) {
      double px = REC_Particle_px->at(i) * REC_Particle_px->at(i);
      double py = REC_Particle_py->at(i) * REC_Particle_py->at(i);
      double pz = REC_Particle_pz->at(i) * REC_Particle_pz->at(i);

      P = TMath::Sqrt(px + py + pz);
      if (P > 0.05) {
        momentum->Fill(P);
      }

      if (REC_Particle_beta->at(i) > 0. && P > 0.0) {
        if (REC_Particle_pid->at(i) != 0) {
          mom_vs_beta_except_e->Fill(P, REC_Particle_beta->at(i));
        }
        if (REC_Particle_charge->at(i) > 0) {
          mom_vs_beta_pos->Fill(P, REC_Particle_beta->at(i));
        } else if (REC_Particle_charge->at(i) < 0) {
          mom_vs_beta_neg->Fill(P, REC_Particle_beta->at(i));
        }
        //  if (REC_Particle_pid->at(i) == 0) {
        //  mom_vs_beta_0->Fill(P, REC_Particle_beta->at(i));
        //}
        if (REC_Particle_pid->at(i) == 2212) {
          mom_vs_beta_proton->Fill(P, REC_Particle_beta->at(i));
        } else if (REC_Particle_pid->at(i) == -2212 &&
                   REC_Particle_charge->at(i) < 0) {
          mom_vs_beta_proton->Fill(P, REC_Particle_beta->at(i));
        } else if (REC_Particle_pid->at(i) == 211) {
          mom_vs_beta_pion->Fill(P, REC_Particle_beta->at(i));
        } else if (REC_Particle_pid->at(i) == -211 &&
                   REC_Particle_charge->at(i) < 0) {
          mom_vs_beta_pion_mi->Fill(P, REC_Particle_beta->at(i));
        } else if (REC_Particle_pid->at(i) == 11) {
          mom_vs_beta_electron->Fill(P, REC_Particle_beta->at(i));
        }
      }

      total++;
    }

    //  }

    /* for index in length(Cal_pindex){
      if(Cal_pindex->at(index) == 0){
  }      energy += Cal_energy->at(index)
      }
      }
      sf = energy / P  */
    // Scintillation bhaneko chamkane ho.
    for (int i = 0; i < REC_Scintillator_time->size(); i++) {
      if (REC_Scintillator_time->size() == 0)
        continue;

      if (REC_Scintillator_detector->at(i) == 12) {
        int j = REC_Scintillator_pindex->at(i);
        if (REC_Particle_pid->at(0) == 11 && REC_Particle_charge->at(0) < 0) {

          //  if (REC_Particle_pid->at(2) == 211)

          //  {
          if (P > 0.0) {

            double electron_vertex =
                vertex_time(REC_Scintillator_time->at(0),
                            REC_Scintillator_path->at(0), 1.0);

            double px = REC_Particle_px->at(j) * REC_Particle_px->at(j);
            double py = REC_Particle_py->at(j) * REC_Particle_py->at(j);
            double pz = REC_Particle_pz->at(j) * REC_Particle_pz->at(j);
            P = TMath::Sqrt(px + py + pz);
            // do for all postitive and all negative particles seperately and
            // then
            // you will see..

            if (REC_Particle_charge->at(j) < 0 && j == 0) {
              double dt_electron_vertex = delta_t(electron_vertex, MASS_E, P,
                                                  REC_Scintillator_time->at(i),
                                                  REC_Scintillator_path->at(i));
              delta_t_electron_vertex_withoutID_vs_P_hist->Fill(
                  P, dt_electron_vertex);
            }
            if (REC_Particle_charge->at(j) < 0 && j == 0 &&
                REC_Particle_pid->at(j) == 11) {
              double dt_electron_vertex_withid = delta_t(
                  electron_vertex, MASS_E, P, REC_Scintillator_time->at(i),
                  REC_Scintillator_path->at(i));
              delta_t_electron_vertex_withID_vs_P_hist->Fill(
                  P, dt_electron_vertex_withid);
            }
            if (REC_Particle_charge->at(j) < 0 && j == 0 &&
                REC_Particle_pid->at(j) != 11) {
              double dt_electron_vertex_anti_cut = delta_t(
                  electron_vertex, MASS_E, P, REC_Scintillator_time->at(i),
                  REC_Scintillator_path->at(i));
              delta_t_electron_vertex_antiID_vs_P_hist->Fill(
                  P, dt_electron_vertex_anti_cut);
            }
            if (REC_Particle_charge->at(j) < 0) { //&&
              //  REC_Particle_pid->at(i) ==
              //    11) { // || REC_Particle_pid->at(j) == -11) {
              double dt_electron = delta_t(electron_vertex, MASS_E, P,
                                           REC_Scintillator_time->at(i),
                                           REC_Scintillator_path->at(i));
              delta_t_electron_vs_P_hist->Fill(P, dt_electron);
            }
            if (REC_Particle_charge->at(j) < 0 &&
                REC_Particle_pid->at(j) != 11) {
              double dt_electron = delta_t(electron_vertex, MASS_E, P,
                                           REC_Scintillator_time->at(i),
                                           REC_Scintillator_path->at(i));
              delta_t_anti_electron_vs_P_hist->Fill(P, dt_electron);
            }
            if (REC_Particle_charge->at(j) < 0 &&
                REC_Particle_pid->at(j) == 11) {
              double dt_electron = delta_t(electron_vertex, MASS_E, P,
                                           REC_Scintillator_time->at(i),
                                           REC_Scintillator_path->at(i));
              delta_t_electron_vs_P_Withid_hist->Fill(P, dt_electron);
            }

            if (REC_Particle_charge->at(j) >
                0) { // REC_Particle_pid->at(j) == 2212) {

              double dt_proton = delta_t(electron_vertex, MASS_P, P,
                                         REC_Scintillator_time->at(i),
                                         REC_Scintillator_path->at(i));

              // if (-10.0 <= dt_proton && dt_proton <= 10.0) {
              delta_t_proton_vs_P_hist->Fill(P, dt_proton);
              //}
            }
            if (REC_Particle_charge->at(j) > 0 &&
                REC_Particle_pid->at(j) == 2212) {

              double dt_proton = delta_t(electron_vertex, MASS_P, P,
                                         REC_Scintillator_time->at(i),
                                         REC_Scintillator_path->at(i));

              //  if (-10.0 <= dt_proton && dt_proton <= 10.0) {
              delta_t_proton_vs_P_Withid_hist->Fill(P, dt_proton);
            }
            if (REC_Particle_charge->at(j) > 0 &&
                REC_Particle_pid->at(j) != 2212) {

              double dt_proton = delta_t(electron_vertex, MASS_P, P,
                                         REC_Scintillator_time->at(i),
                                         REC_Scintillator_path->at(i));
              delta_t_anti_proton_vs_P_hist->Fill(P, dt_proton);
              if (P > 0) {
                proton_P_hist->Fill(P);
              }
            }
            if (REC_Particle_charge->at(j) >
                0) { // REC_Particle_pid->at(j) == 211) {
              double dt_pip = delta_t(electron_vertex, MASS_PIP, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_pip_vs_P_hist->Fill(P, dt_pip);
            }
            if (REC_Particle_charge->at(j) > 0 &&
                REC_Particle_pid->at(j) == 211) {
              double dt_pip = delta_t(electron_vertex, MASS_PIP, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_pip_vs_P_Withid_hist->Fill(P, dt_pip);
            }
            if (REC_Particle_charge->at(j) > 0 &&
                REC_Particle_pid->at(j) != 211) {
              double dt_pip = delta_t(electron_vertex, MASS_PIP, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_anti_pip_vs_P_hist->Fill(P, dt_pip);
            }
            if (REC_Particle_charge->at(j) >
                0) { // REC_Particle_pid->at(j) == 321) {
              double dt_kap = delta_t(electron_vertex, MASS_KP, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_kap_vs_P_hist->Fill(P, dt_kap);
            }

            if (REC_Particle_charge->at(j) > 0 &&
                REC_Particle_pid->at(j) == 321) {
              double dt_kap = delta_t(electron_vertex, MASS_KP, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_kap_vs_P_Withid_hist->Fill(P, dt_kap);
            }
            if (REC_Particle_charge->at(j) < 0 &&
                REC_Particle_pid->at(j) != 11) {
              double dt_pim = delta_t(electron_vertex, MASS_PIM, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_pim_vs_P_hist->Fill(P, dt_pim);
            }
            if (REC_Particle_charge->at(j) < 0 &&
                REC_Particle_pid->at(j) == -211) {
              double dt_pim = delta_t(electron_vertex, MASS_PIM, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_pim_vs_P_Withid_hist->Fill(P, dt_pim);
            }
            if (REC_Particle_charge->at(j) ==
                0) { // REC_Particle_pid->at(j) == 111) {
              double dt_pi0 = delta_t(electron_vertex, MASS_PI0, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_pi0_vs_P_hist->Fill(P, dt_pi0);
            }
            if (REC_Particle_charge->at(j) == 0 &&
                REC_Particle_pid->at(j) == 111) {
              double dt_pi0 = delta_t(electron_vertex, MASS_PI0, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_pi0_vs_P_hist->Fill(P, dt_pi0);
            }
            if (REC_Particle_charge->at(j) < 0 &&
                REC_Particle_pid->at(j) != 11) {
              double dt_kam = delta_t(electron_vertex, MASS_KM, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_kam_vs_P_hist->Fill(P, dt_kam);
            }
            if (REC_Particle_charge->at(j) < 0 &&
                REC_Particle_pid->at(j) == -321) {
              double dt_kam = delta_t(electron_vertex, MASS_KM, P,
                                      REC_Scintillator_time->at(i),
                                      REC_Scintillator_path->at(i));
              delta_t_kam_vs_P_Withid_hist->Fill(P, dt_kam);
            }
            if (REC_Particle_charge->at(j) ==
                0) { // REC_Particle_pid->at(j) == 22) {
              double dt_g = delta_t(electron_vertex, MASS_G, P,
                                    REC_Scintillator_time->at(i),
                                    REC_Scintillator_path->at(i));
              delta_t_g_vs_P_hist->Fill(P, dt_g);
            }
            if (REC_Particle_charge->at(j) == 0 &&
                REC_Particle_pid->at(j) == 22) {
              double dt_g = delta_t(electron_vertex, MASS_G, P,
                                    REC_Scintillator_time->at(i),
                                    REC_Scintillator_path->at(i));
              delta_t_g_vs_P_Withid_hist->Fill(P, dt_g);
            }

            // delta_t_pip_vs_P_hist->Fill(P, dt_pip);
            // delta_t_kap_vs_P_hist->Fill(P, dt_kap);
            //  delta_t_proton_vs_P_hist->Fill(P, dt_proton);
            //    delta_t_electron_vs_P_hist->Fill(P, dt_electron);
            //  delta_t_pi0_vs_P_hist->Fill(P, dt_pi0);
            //  delta_t_kam_vs_P_hist->Fill(P, dt_kam);
            //  delta_t_pim_vs_P_hist->Fill(P, dt_pim);
            //  delta_t_g_vs_P_hist->Fill(P, dt_g);
            //    delta_t_omega_vs_P_hist->Fill(P, dt_omega);
          }
        }
      }
    }
    if (REC_Particle_pid->at(0) == 11) {
      // Setup scattered electron 4 vector
      TVector3 e_mu_prime_3;
      TLorentzVector e_mu_prime;
      TLorentzVector e_mu(0.0, 0.0, 2.2, 2.2);
      e_mu_prime_3.SetXYZ(REC_Particle_px->at(0), REC_Particle_py->at(0),
                          REC_Particle_pz->at(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      double W = W_calc(e_mu, e_mu_prime);
      double Q2 = Q2_calc(e_mu, e_mu_prime);

      W_hist->Fill(W);
      //  if (Q2 > 0.5 && Q2 < 5.5) {
      Q2_hist->Fill(Q2);
      //  }
      if (Q2 > 0.0) {
        W_vs_q2->Fill(W, Q2);
      }
    }
    //  }

    for (int i = 0; i < REC_Traj_pindex->size(); i++) {
      if (REC_Traj_pindex->size() == 0)
        continue;

      // std::cout << "my det dc " << REC_Traj_x->at(index);
      int j = REC_Traj_pindex->at(i);
      if ( //(REC_Particle_pid->at(j) == 321) || // && //j == 0)
          (REC_Particle_pid->at(j) == -2212) //||
          //  (REC_Particle_pid->at(j) == 211)
          //(REC_Particle_pid->at(j) == 2212)
          && (REC_Particle_charge->at(j) < 0)) {
        if (REC_Traj_detId->at(i) == 24) {
          double angle = 60;
          double height = 38;
          for (int index = 0; index < REC_Track_pindex->size(); index++) {
            int sec = REC_Track_sector->at(index) - 1;
            double x1_rot = REC_Traj_y->at(i) * sin(sec * 60.0 * Pival / 180) +
                            REC_Traj_x->at(i) * cos(sec * 60.0 * Pival / 180);
            double y1_rot = REC_Traj_y->at(i) * cos(sec * 60.0 * Pival / 180) -
                            REC_Traj_x->at(i) * sin(sec * 60.0 * Pival / 180);
            double slope = 1 / tan(0.5 * angle * Pival / 180);
            double left = (height - slope * y1_rot);
            double right = (height + slope * y1_rot);
            double radius2_DCr1 = pow(46, 2) - pow(y1_rot, 2);
            if (x1_rot > left && x1_rot > right &&
                pow(x1_rot, 2) > radius2_DCr1) {
              dc_Fiducial_cut_hist->Fill(REC_Traj_x->at(i), REC_Traj_y->at(i));
            }
          }
        }
      }
    }
    /*  std::vector<float> *final_x_coord = new std::vector<float>();
      std::vector<float> *final_y_coord = new std::vector<float>();

      // you tringular cut ho, stafan ko rotation saga farak chha.

      for (int i = 0; i < REC_Calorimeter_pindex->size(); i++) {
        if (REC_Calorimeter_pindex->size() == 0)
          continue;

        int j = REC_Calorimeter_pindex->at(i);
        if (REC_Particle_pid->at(j) == 11 && j == 0)

          if (REC_Calorimeter_sector->at(i) == 1 &&
              REC_Calorimeter_x->at(i) < 350) {
            // printf(" Here I am\n");
            final_x_coord->push_back(REC_Calorimeter_x->at(i));
            final_y_coord->push_back(REC_Calorimeter_y->at(i));
          }

          else if (REC_Calorimeter_sector->at(i) == 2 &&
                   REC_Calorimeter_x->at(i) > 0 &&
                   REC_Calorimeter_y->at(i) > 0.577 * REC_Calorimeter_x->at(i)
      &&
                   REC_Calorimeter_y->at(i) <
                       362 - (0.577 * REC_Calorimeter_x->at(i)))

          {
            final_x_coord->push_back(REC_Calorimeter_x->at(i));
            final_y_coord->push_back(REC_Calorimeter_y->at(i));

          } else if (REC_Calorimeter_sector->at(i) == 3 &&
                     REC_Calorimeter_x->at(i) < -0.5 &&
                     REC_Calorimeter_y->at(i) >
                         -0.577 * REC_Calorimeter_x->at(i) &&
                     REC_Calorimeter_y->at(i) <
                         362 + (0.577 * REC_Calorimeter_x->at(i)))

          {
            final_x_coord->push_back(REC_Calorimeter_x->at(i));
            final_y_coord->push_back(REC_Calorimeter_y->at(i));
          } else if (REC_Calorimeter_sector->at(i) == 4 &&
                     REC_Calorimeter_x->at(i) > -362 &&
                     REC_Calorimeter_y->at(i) <=
                         -0.53 * (REC_Calorimeter_x->at(i) - 44))

          {
            final_x_coord->push_back(REC_Calorimeter_x->at(i));
            final_y_coord->push_back(REC_Calorimeter_y->at(i));
          } else if (REC_Calorimeter_sector->at(i) == 5 &&
                     REC_Calorimeter_x->at(i) < 0 &&
                     REC_Calorimeter_y->at(i) <
                         0.577 * REC_Calorimeter_x->at(i) &&
                     REC_Calorimeter_y->at(i) >
                         -362 - (0.577 * REC_Calorimeter_x->at(i)))

          {
            final_x_coord->push_back(REC_Calorimeter_x->at(i));
            final_y_coord->push_back(REC_Calorimeter_y->at(i));
          } else if (REC_Calorimeter_sector->at(i) == 6 &&
                     REC_Calorimeter_x->at(i) > 0 &&
                     REC_Calorimeter_y->at(i) <
                         -0.577 * REC_Calorimeter_x->at(i) &&
                     REC_Calorimeter_y->at(i) >
                         -362 + (0.577 * REC_Calorimeter_x->at(i)))

          {
            final_x_coord->push_back(REC_Calorimeter_x->at(i));
            final_y_coord->push_back(REC_Calorimeter_y->at(i));
          }

        for (int j = 0; j < final_x_coord->size(); j++) {
          pcal_Fiducial_cut_hist->Fill( // x_PCAL, y_PCAL);
              final_x_coord->at(j), final_y_coord->at(j));
          //  final_x_coord->clear();
          //  final_y_coord->clear();
        }
      }*/
    double energy = 0;
    double PCAL = 0;
    double ECAL = 0;
    double ECin = 0;
    double ECout = 0;
    TLorentzVector e_mu_prime;
    TVector3 e_mu_prime_3;

    for (int i = 0; i < REC_Calorimeter_pindex->size(); i++) {
      if (REC_Calorimeter_pindex->size() == 0)
        continue;
      int j = REC_Calorimeter_pindex->at(i);
      if (REC_Particle_pid->at(j) == 11 && j == 0) {

        // bool EC_hit_position_fiducial_cut(int j)

        int sec_PCAL = REC_Calorimeter_sector->at(i) - 1; //
        // part_Cal_PCAL_sector[j] - 1;
        double x_PCAL = REC_Calorimeter_x->at(i);
        //  part_Cal_PCAL_x[j];
        double y_PCAL = REC_Calorimeter_y->at(i);
        //  part_Cal_PCAL_y[j];
        double x_PCAL_rot = y_PCAL * sin(sec_PCAL * 60.0 * Pival / 180) +
                            x_PCAL * cos(sec_PCAL * 60.0 * Pival / 180);
        double y_PCAL_rot = y_PCAL * cos(sec_PCAL * 60.0 * Pival / 180) -
                            x_PCAL * sin(sec_PCAL * 60.0 * Pival / 180);
        double angle_PCAL = 60;
        double height_PCAL = 45;
        double slope_PCAL =
            1 / tan(0.5 * angle_PCAL * Pival / 180); // about =1.73
        double left_PCAL =
            (height_PCAL - slope_PCAL * y_PCAL_rot); // about -214 for
                                                     //    150

        double right_PCAL =
            (height_PCAL + slope_PCAL * y_PCAL_rot); // about 304 for
                                                     //    150
        double radius2_PCAL = pow(height_PCAL + 6, 2) - pow(y_PCAL_rot, 2);
        // 51*51 - 150*150 = 19,899
        if (x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL &&
            pow(x_PCAL_rot, 2) > radius2_PCAL && x_PCAL_rot < 362)

        {
          pcal_Fiducial_cut_hist->Fill(x_PCAL, y_PCAL);
        }
        //}
        //  }
        //  return true;
        //  else
        //  return false;

        e_mu_prime_3.SetXYZ(REC_Particle_px->at(j), REC_Particle_py->at(j),
                            REC_Particle_pz->at(j));
        P = e_mu_prime_3.Mag();
        e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

        //  for (int i = 0; i < REC_Calorimeter_pindex->size(); i++) {
        //  if (REC_Calorimeter_sector->at(i) == 1) {
        //    if (REC_Calorimeter_index->at(i) == 0) {
        //  if (REC_Calorimeter_sector->at(i) == 4) {
        energy = energy + REC_Calorimeter_energy->at(i);

        if (REC_Calorimeter_layer->at(i) == 1) {
          PCAL = REC_Calorimeter_energy->at(i);
        }
        if (REC_Calorimeter_layer->at(i) == 4) {
          ECin = REC_Calorimeter_energy->at(i);
        }
        if (REC_Calorimeter_layer->at(i) == 7) {
          ECout = REC_Calorimeter_energy->at(i);
        }

        //}
      }

      { total++; }
    }
    SF = energy / e_mu_prime.P();
    //  if (SF > 0.05 && SF < 0.6) {
    if (energy != 0) {
      SF_vs_P_hist->Fill(e_mu_prime.P(), SF);
    }
    ECAL = ECin + ECout;
    if (PCAL != 0) {
      PCAL_hist->Fill(PCAL);
    }
    if (ECAL != 0) {
      ECAL_hist->Fill(ECAL);
    }
    if (ECAL != 0) {
      if (PCAL > 0.01)
        ECAL_vs_PCAL_hist->Fill(PCAL, ECAL);
    }
    if (ECout != 0) {
      if (ECin > 0.01)
        ECout_vs_ECin_hist->Fill(ECin, ECout);
    }

    pcal_SF = PCAL / e_mu_prime.P();
    if (pcal_SF != 0) {
      pcal_SF_vs_P_hist->Fill(e_mu_prime.P(), pcal_SF);
    }
    //  }
    ecal_SF = ECAL / e_mu_prime.P();
    //  if (ecal_SF > 0.05 && ecal_SF < 0.8) {
    if (ecal_SF != 0) {
      ecal_SF_vs_P_hist->Fill(e_mu_prime.P(), ecal_SF);
    }
  }
  //}

  out->cd();
  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  // wvsq2->GetOption("COLZ");
  momentum->Write();
  W_hist->Write();
  Q2_hist->Write();
  W_vs_q2->Write();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  //  mom_vs_beta->GetOption("COLZ");
  mom_vs_beta_except_e->Write();
  mom_vs_beta_pos->Write();
  mom_vs_beta_neg->Write();
  //  mom_vs_beta_neg->Fit("gaus");
  mom_vs_beta_proton->Write();
  mom_vs_beta_anti_proton->Write();
  mom_vs_beta_pion->Write();
  mom_vs_beta_pion_mi->Write();
  mom_vs_beta_electron->Write();
  //  mom_vs_beta_0->Write();

  TDirectory *SF_VS_P = out->mkdir("SF_vs_P");
  SF_VS_P->cd();
  //  SF_VS_P->GetOption("COLZ");
  SF_vs_P_hist->Write();
  PCAL_hist->Write();
  ECAL_hist->Write();
  pcal_SF_vs_P_hist->Write();
  ecal_SF_vs_P_hist->Write();
  ECAL_vs_PCAL_hist->Write();
  ECout_vs_ECin_hist->Write();
  pcal_Fiducial_cut_hist->Write();
  Z_Vertex_position_hist->Write();
  dc_Fiducial_cut_hist->Write();

  TDirectory *delta_t_vs_P = out->mkdir("delta_t_vs_P");
  delta_t_vs_P->cd();
  // delta_t_vs_P->GetOption("COLZ");
  delta_t_electron_vertex_withoutID_vs_P_hist->Write();
  delta_t_electron_vertex_withID_vs_P_hist->Write();
  delta_t_electron_vertex_antiID_vs_P_hist->Write();
  delta_t_electron_vs_P_hist->Write();
  delta_t_anti_electron_vs_P_hist->Write();
  delta_t_electron_vs_P_Withid_hist->Write();
  delta_t_proton_vs_P_hist->Write();
  proton_P_hist->Write();
  delta_t_anti_proton_vs_P_hist->Write();
  delta_t_proton_vs_P_Withid_hist->Write();
  delta_t_pip_vs_P_hist->Write();
  delta_t_anti_pip_vs_P_hist->Write();
  delta_t_pip_vs_P_Withid_hist->Write();
  delta_t_kap_vs_P_hist->Write();
  delta_t_kap_vs_P_Withid_hist->Write();
  delta_t_pi0_vs_P_hist->Write();
  delta_t_pi0_vs_P_Withid_hist->Write();
  delta_t_pim_vs_P_hist->Write();
  delta_t_pim_vs_P_Withid_hist->Write();
  delta_t_kam_vs_P_hist->Write();
  delta_t_kam_vs_P_Withid_hist->Write();
  delta_t_g_vs_P_hist->Write();
  delta_t_g_vs_P_Withid_hist->Write();
  //  delta_t_omega_vs_P_hist->Write();

  out->Close();
  chain.Reset();
  std::cerr << "\n" << total << std::endl;
}

#endif
