//
// Created by yuehu on 2021/9/8.
//
//     int nband;
//     using BandList = vector<Band>;
//     BandList bands;
//     vector<int> geo_matter_index;
//     vector<double> kappabulk;

#include <sstream>
#include <cmath>
#include "BTEBand/BTEBand.h"
using namespace std;
BTEBand::BTEBand(ifstream &inFile, int Dimension_Material)
{
    string ant="!";
    if (!inFile.is_open())
    {
        cout << "DEBUG: band file not open" << endl;
        exit(1);
    }
    string str;
    string line;
    int Num_Matter=0;
    while(getline(inFile, str))
    {
        if(str.find("MaterialNumber") >= 0 && str.find("MaterialNumber") < str.length() && str.find(ant) == str.npos)
        {
            inFile >> Num_Matter;
        }

        if(str.find("BandNumber") >= 0 && str.find("BandNumber") < str.length() && str.find(ant) == str.npos)
        {
            inFile >> nband;
        }
    }
    if (Num_Matter==0||nband==0)
    {
        cout<<"Error: please provide Number of matter and Number of bands"<<endl;
    }

    char new_line;
    bands.resize(Num_Matter);
    bands.resize(nband);
    for (int i = 0; i < nband; ++i)
    {
        bands[i].group_velocity.resize(Num_Matter);
        bands[i].relaxation_time.resize(Num_Matter);
        bands[i].heat_capacity.resize(Num_Matter);
        bands[i].lattice_ratio.resize(Num_Matter);
        bands[i].heat_ratio.resize(Num_Matter);
    }
    kappabulk.resize(Num_Matter);
    capacitybulk.resize(Num_Matter);
    inFile.close();
    inFile.open("input/PHONON");
    for (int i = 0; i < Num_Matter; ++i)
    {
        while (getline(inFile, line))
        {
            if (line.find("MaterialFile") >= 0 && line.find("MaterialFile") < line.length() && str.find(ant) == str.npos)
            {

                if (line.find("DATABASE") > 0 && line.find("DATABASE") < line.length())
                {
                    string File;
                    inFile>>File;
                    ifstream inputband(File);
                    if (!inputband.is_open())
                    {
                        cout << "Error: DATABASE band file not open" << endl;
                        exit(1);
                    }
                    for (int j = 0; j < nband; ++j)
                    {
                        inputband>> bands[j].group_velocity[i] >> bands[j].relaxation_time[i] >> bands[j].heat_capacity[i] >> bands[j].lattice_ratio[i] >> bands[j].heat_ratio[i];
                    }


                }
                else if (line.find("ShengBTE") > 0 && line.find("ShengBTE") < line.length())
                {
                    string Directory;
                    inFile>>Directory;
                    double ngrid;
                    double lfactor=1.0;
                    double lattvec[3][3];
                    double Temperature;
                    double Volume;

                    string File_control=Directory + "/CONTROL";//?
                    string tempstr;
                    ifstream control(File_control);
                    if (!control.is_open())
                    {
                        cout << "Error: ShengBTE CONTROL file not open" << endl;
                        exit(1);
                    }
                    while (getline(control, tempstr))
                    {
                        if (tempstr.find("ngrid(:)") > 0 && tempstr.find("ngrid(:)") < tempstr.length())
                        {
                            string vv1 = tempstr.substr(tempstr.find("=")+1);
                            string vv2;
                            int pos1=0;
                            int pos2=0;
                            for (int si=0; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            ngrid = stod(vv1.substr(pos1, pos2-pos1));
                        }
                        if (tempstr.find("lfactor") > 0 && tempstr.find("lfactor") < tempstr.length())
                        {
                            string vv1 = tempstr.substr(tempstr.find("=")+1);
                            string vv2;
                            int pos1=0;
                            for (int si=0; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            lfactor = stod(vv1.substr(pos1, vv1.find(",")));
                        }
                        if (tempstr.find("lattvec(:,1)") > 0 && tempstr.find("lattvec(:,1)") < tempstr.length())
                        {
                            string vv1 = tempstr.substr(tempstr.find("=")+1);
                            string vv2;
                            int pos1=0;
                            int pos2=0;
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            lattvec[0][0] = stod(vv1.substr(pos1, pos2-pos1));

                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            lattvec[0][1] = stod(vv1.substr(pos1, pos2-pos1));

                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            lattvec[0][2] = stod(vv1.substr(pos1));
                        }
                        if (tempstr.find("lattvec(:,2)") > 0 && tempstr.find("lattvec(:,2)") < tempstr.length())
                        {
                            string vv1 = tempstr.substr(tempstr.find("=")+1);
                            string vv2;
                            int pos1=0;
                            int pos2=0;
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            lattvec[1][0] = stod(vv1.substr(pos1, pos2-pos1));

                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            lattvec[1][1] = stod(vv1.substr(pos1, pos2-pos1));

                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            lattvec[1][2] = stod(vv1.substr(pos1));
                        }
                        if (tempstr.find("lattvec(:,3)") > 0 && tempstr.find("lattvec(:,3)") < tempstr.length())
                        {
                            string vv1 = tempstr.substr(tempstr.find("=")+1);
                            string vv2;
                            int pos1=0;
                            int pos2=0;
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            lattvec[2][0] = stod(vv1.substr(pos1, pos2-pos1));

                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            lattvec[2][1] = stod(vv1.substr(pos1, pos2-pos1));

                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            lattvec[2][2] = stod(vv1.substr(pos1));
                        }
                        if (tempstr.find("T") > 0 && tempstr.find("T") < tempstr.length())
                        {
                            string vv1 = tempstr.substr(tempstr.find("=")+1);
                            string vv2;
                            int pos1=0;
                            for (int si=0; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            //Temperature = stod(vv1.substr(pos1));//unit:K
                            if (vv1.find("TRUE")==str.npos && vv1.find("True")==str.npos)
                            {
                                Temperature = stod(vv1.substr(pos1));//unit:K
                            }
                        }
                    }
                    Volume = (lattvec[0][0]*lattvec[1][1]*lattvec[2][2] + lattvec[1][0]*lattvec[2][1]*lattvec[0][2] + lattvec[0][1]*lattvec[1][2]*lattvec[2][0] - lattvec[0][2]*lattvec[1][1]*lattvec[2][0] - lattvec[0][1]*lattvec[0][1]*lattvec[2][2] - lattvec[0][0]*lattvec[1][2]*lattvec[2][1])*lfactor*lfactor*lfactor;//unit:?

                    vector<vector<double>> omegas;
                    vector<double> omegas1;
                    string File_omega=Directory+"/BTE.omega";//?
                    ifstream omega(File_omega);
                    if (!omega.is_open())
                    {
                        cout << "Error: ShengBTE BTE.omega file not open" << endl;
                        exit(1);
                    }
                    int ioo=0;
                    int joo;
                    while (getline(omega, tempstr))
                    {
                        if (tempstr.length() <= 1)
                        {
                            break;
                        }
                        string vv1 = tempstr;
                        string vv2;
                        int pos1=0;
                        int pos2=0;
                        joo=0;
                        while (pos2>=pos1)
                        {
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            if(pos2>=pos1){
                                double value=stod(vv1.substr(pos1, pos2-pos1));//unit:rad/ps
                                omegas1.push_back(value);
                            }else{
                                double value=stod(vv1.substr(pos1, vv1.size()-pos1));
                                omegas1.push_back(value);
                            }
                            joo++;
                        }
                        omegas.push_back(omegas1);
                        omegas1.clear();
                        ioo++;
                    }

                    vector<vector<double>> qpoint;
                    vector<double> qpoint1;
                    string File_qpoints=Directory+"/BTE.qpoints_full";//?
                    ifstream qpoints(File_qpoints);
                    if (!qpoints.is_open())
                    {
                        cout << "DEBUG: ShengBTE BTE.qpoints_full file not open" << endl;
                        exit(1);
                    }
                    int iqq=0;
                    int jqq;
                    while (getline(qpoints, tempstr))
                    {
                        if (tempstr.length() <= 1)
                        {
                            break;
                        }
                        string vv1 = tempstr;
                        string vv2;
                        int pos1=0;
                        int pos2=0;
                        jqq=0;
                        while (pos2>=pos1)
                        {
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            if(pos2>=pos1){
                                double value=stod(vv1.substr(pos1, pos2-pos1));
                                qpoint1.push_back(value);
                            }else{
                                double value=stod(vv1.substr(pos1, vv1.size()-pos1));
                                qpoint1.push_back(value);
                            }
                            jqq++;
                        }
                        qpoint.push_back(qpoint1);
                        qpoint1.clear();
                        iqq++;
                    }

                    vector<vector<double>> w_final;
                    vector<double> w_final1;
                    string File_wfinal=Directory+"/BTE.w_final";//?
                    ifstream wfinal(File_wfinal);
                    if (!wfinal.is_open())
                    {
                        cout << "Error: ShengBTE BTE.w_final file not open" << endl;
                        exit(1);
                    }
                    int iww=0;
                    int jww;
                    while (getline(wfinal, tempstr))
                    {
                        if (tempstr.length() <= 1)
                        {
                            break;
                        }
                        string vv1 = tempstr;
                        string vv2;
                        int pos1=0;
                        int pos2=0;
                        jww=0;
                        while (pos2>=pos1)
                        {
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            if(pos2>=pos1){
                                double value=stod(vv1.substr(pos1, pos2-pos1));//unit:1/ps
                                w_final1.push_back(value);
                            }else{
                                double value=stod(vv1.substr(pos1, vv1.size()-pos1));//unit:1/ps
                                w_final1.push_back(value);
                            }
                            jww++;
                        }
                        w_final.push_back(w_final1);
                        w_final1.clear();
                        iww++;
                    }

                    vector<vector<double>> v_full;
                    vector<double> v_full1;
                    string File_vfull=Directory+"/BTE.v";//?
                    ifstream vfull(File_vfull);
                    if (!vfull.is_open())
                    {
                        cout << "DEBUG: ShengBTE BTE.v file not open" << endl;
                        exit(1);
                    }
                    int ivv=0;
                    int jvv;
                    while (getline(vfull, tempstr))
                    {
                        if (tempstr.length() <= 1)
                        {
                            break;
                        }
                        string vv1 = tempstr;
                        string vv2;
                        int pos1=0;
                        int pos2=0;
                        jvv=0;
                        while (pos2>=pos1)
                        {
                            for (int si=pos2; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2!=" "){
                                    pos1 = si;
                                    break;
                                }
                            }
                            for (int si=pos1; si<vv1.size(); si++){
                                vv2 = vv1.substr(si,1);
                                if (vv2==" "){
                                    pos2 = si;
                                    break;
                                }
                            }
                            if(pos2>=pos1){
                                double value=stod(vv1.substr(pos1, pos2-pos1));//unit:km/s
                                v_full1.push_back(value);
                            }else{
                                double value=stod(vv1.substr(pos1, vv1.size()-pos1));//unit:km/s
                                v_full1.push_back(value);
                            }
                            jvv++;
                        }
                        v_full.push_back(v_full1);
                        v_full1.clear();
                        ivv++;
                    }

                    vector<vector<double>> omega_all, tau_all, vel_all, capacity_all;
                    vector<double> temp1, temp2, temp3, temp4;
                    for (int i1=0; i1<iqq; i1++)
                    {
                        for (int i2=0; i2<joo; i2++)
                        {
                            temp1.push_back(omegas[qpoint[i1][1]-1][i2] * 1E12); //unit: rad/s
                            temp2.push_back(1/w_final[i2*ioo+qpoint[i1][1]-1][1] * 1E-12);//unit: s
                            //temp2.push_back(w_final[i2*ioo+qpoint[i1][1]-1][1]);//unit: 1/ps
                            temp3.push_back(pow(pow(v_full[i2*ioo+qpoint[i1][1]-1][0]*1000,2) + pow(v_full[i2*ioo+qpoint[i1][1]-1][1]*1000,2) + pow(v_full[i2*ioo+qpoint[i1][1]-1][2]*1000,2),0.5));//unit: m/s
                        }
                        omega_all.push_back(temp1);//nqpoints*nbands
                        tau_all.push_back(temp2);
                        vel_all.push_back(temp3);
                        temp1.clear();
                        temp2.clear();
                        temp3.clear();
                    }

                    double h=6.62620E-34;//unit:Js
                    double kb=1.38062E-23;//unit:J/K
                    double hbar=1.05459E-34;//unit:Js
                    Volume=abs(Volume)*1E-27;//unit:nm^-3
                    for (int i1=0; i1<iqq; i1++)
                    {
                        for (int i2=0; i2<joo; i2++)
                        {
                            double con = hbar * omega_all[i1][i2] / kb / Temperature;
                            double value = 1 / Volume * kb * con * con * exp(con) / pow(exp(con)-1,2);
                            temp4.push_back(value);
                        }
                        capacity_all.push_back(temp4);
                        temp4.clear();
                    }

                    vector<double> tau_mode, vel_mode, capacity_mode, kappa_mode, frequency_mode;// nmode*1
                    int num_mode=0;
                    int num_nan=0;
                    for (int i2=0; i2<joo; i2++)
                    {
                        for (int i1=0; i1<iqq; i1++)
                        {
                            double value = capacity_all[i1][i2] / iqq * vel_all[i1][i2] * vel_all[i1][i2] * tau_all[i1][i2];
                            if (isnan(value)){
                                num_nan++;
                            }else if(value>=0){
                                frequency_mode.push_back(omega_all[i1][i2]);
                                capacity_mode.push_back(capacity_all[i1][i2] / iqq);
                                vel_mode.push_back(vel_all[i1][i2]);
                                tau_mode.push_back(tau_all[i1][i2]);
                                kappa_mode.push_back(value);
                                num_mode++;
                            }
                        }
                    }

                    int totalmode=frequency_mode.size();
                    double maxfreq=0;
                    double minfreq=1e20;
                    for (int j = 0; j < totalmode; ++j) {
                        if(frequency_mode[j]>maxfreq)
                        {
                            maxfreq=frequency_mode[j];
                        }
                        if(frequency_mode[j]< minfreq)
                        {
                            minfreq=frequency_mode[j];
                        }
                    }
                    vector<double> freqbin;
                    freqbin.resize(nband+1);
                    freqbin[0]=-1;//
                    for (int j = 1; j < nband+1; ++j) {
                        //freqbin[j]=minfreq+(maxfreq-minfreq)/nband*j;
                        freqbin[j]=freqbin[0]+ceil(totalmode/nband)*j;
                        if (j==nband)
                        {
                            freqbin[j]=totalmode;
                        }
                    }

                    vector<double> sumC,sumCv,sumk;
                    sumC.resize(nband,0);
                    sumCv.resize(nband,0);
                    sumk.resize(nband,0);
                    for (int j = 0; j < totalmode; ++j) {
                        for (int k = 0; k < nband; ++k) {
                            if(j>freqbin[k]&&j<=freqbin[k+1])
                            {
                                sumC[k]+=capacity_mode[j];
                                sumCv[k]+=capacity_mode[j]*vel_mode[j];
                                sumk[k]+=kappa_mode[j];
                            }
                        }
                    }
                    double sumL=0;
                    double sumheat=0;
                    for (int j = 0; j < nband; ++j) {
                        bands[j].group_velocity[i]=sumCv[j]/sumC[j];
                        bands[j].heat_capacity[i]=sumC[j];
                        bands[j].relaxation_time[i]=sumk[j]/bands[j].heat_capacity[i]/bands[j].group_velocity[i]/bands[j].group_velocity[i];
                        sumL+= bands[j].heat_capacity[i]/ bands[j].relaxation_time[i];
                        sumheat+=bands[j].heat_capacity[i];
                    }

                    for (int j = 0; j < nband; ++j) {
                        bands[j].lattice_ratio[i]= bands[j].heat_capacity[i]/ bands[j].relaxation_time[i]/sumL;
                        bands[j].heat_ratio[i]=bands[j].heat_capacity[i]/sumheat;
                    }
                    double allkappa=0;
                    double allheat=0;
                    for (int j = 0; j < nband; ++j) {
                        allkappa+=bands[j].heat_capacity[i]*bands[j].group_velocity[i]*bands[j].group_velocity[i]*bands[j].relaxation_time[i];
                        allheat+=bands[j].heat_capacity[i];
                    }
                    //cout<<endl;

                    for (int j = 0; j < nband; ++j)
                    {
                        //cout<<"Material #"<< i <<" band #"<< j << " " <<bands[j].group_velocity[i]<<" "<<bands[j].relaxation_time[i]<<" "<<bands[j].heat_capacity[i]<<" "<<bands[j].heat_capacity[i]*bands[j].relaxation_time[i]*bands[j].group_velocity[i]*bands[j].group_velocity[i]<<endl;
                    }
                }
                else if (line.find("ALAMODE") > 0 && line.find("ALAMODE") < line.length())
                {
                    string File;
                    inFile>>File;
                    ifstream inputband(File);
                    if (!inputband.is_open())
                    {
                        cout << "Error: ALAMODE band file not open" << endl;
                        exit(1);
                    }
                    string strband;
                    while (getline(inputband, strband))
                    {

                        if(strband.find("Phonon lifetime")>0 && strband.find("Phonon lifetime")<line.length())
                        {
                            string a1,a2;
                            int startq,endq;
                            int startbranch,endbranch;
                            inputband>>a1>>a2>>a1>>startq>>endq;
                            inputband>>a1>>a2>>a1>>startbranch>>endbranch;
                            getline(inputband, strband);
                            getline(inputband, strband);
                            vector<double> frequency,tau,velocity,mfp,multici,kappa_x,kappa_y,kappa_z,kappa,capacity;
                            int totalmode=(endbranch-startbranch+1)*(endq-startq+1);
                            frequency.resize(totalmode);
                            velocity.resize(totalmode);
                            tau.resize(totalmode);
                            mfp.resize(totalmode);
                            multici.resize(totalmode);
                            kappa_x.resize(totalmode);
                            kappa_y.resize(totalmode);
                            kappa_z.resize(totalmode);
                            kappa.resize(totalmode);
                            capacity.resize(totalmode);
                            for (int j = 0; j < totalmode; ++j) {
                                int indexq,indexb;
                                inputband>>indexq>>indexb>>frequency[j]>>tau[j]>>velocity[j]>>mfp[j]>>multici[j]>>kappa_x[j]>>kappa_y[j]>>kappa_y[j]>>kappa_y[j]>>kappa_y[j]>>kappa_z[j]>>kappa_z[j]>>kappa_z[j]>>kappa_z[j];
                            }

                            for (int j = 0; j < totalmode; ++j) {
                                tau[j]=tau[j]*1e-12;
                                kappa[j]=(kappa_x[j]+kappa_y[j]+kappa_z[j])*multici[j];
                                capacity[j]=kappa[j]/velocity[j]/velocity[j]/tau[j];
                                if(tau[j]==0||velocity[j]==0)
                                {
                                    capacity[j]=0;
                                }
                            }
                            double maxfreq=0;
                            double minfreq=1e20;
                            for (int j = 0; j < totalmode; ++j) {
                                if(frequency[j]>maxfreq)
                                {
                                    maxfreq=frequency[j];
                                }
                                if(frequency[j]< minfreq)
                                {
                                    minfreq= frequency[j];
                                }
                            }
                            vector<double> freqbin;
                            freqbin.resize(nband+1);
                            freqbin[0]=-1;
                            for (int j = 1; j < nband+1; ++j) {
                                freqbin[j]=freqbin[0]+ceil(totalmode/nband)*j;
                                if (j==nband)
                                {
                                    freqbin[j]=totalmode;
                                }
                            }
                            vector<double> sumC,sumCv,sumk;
                            sumC.resize(nband,0);
                            sumCv.resize(nband,0);
                            sumk.resize(nband,0);
                            for (int j = 0; j < totalmode; ++j) {
                                for (int k = 0; k < nband; ++k) {
                                    if(j>freqbin[k]&&j<=freqbin[k+1])
                                    {
                                        sumC[k]+=capacity[j];
                                        sumCv[k]+=capacity[j]*velocity[j];
                                        sumk[k]+=kappa[j];
                                    }
                                }
                            }
                            double sumL=0;
                            double sumheat=0;
                            for (int j = 0; j < nband; ++j) {
                                bands[j].group_velocity[i]=sumCv[j]/sumC[j];
                                bands[j].heat_capacity[i]=sumC[j];
                                bands[j].relaxation_time[i]=sumk[j]/bands[j].heat_capacity[i]/bands[j].group_velocity[i]/bands[j].group_velocity[i];
                                sumL+= bands[j].heat_capacity[i]/ bands[j].relaxation_time[i];
                                sumheat+=bands[j].heat_capacity[i];
                            }

                            for (int j = 0; j < nband; ++j) {
                                bands[j].lattice_ratio[i]= bands[j].heat_capacity[i]/ bands[j].relaxation_time[i]/sumL;
                                bands[j].heat_ratio[i]=bands[j].heat_capacity[i]/sumheat;
                            }
                            double allkappa=0;
                            double allheat=0;
                            for (int j = 0; j < nband; ++j) {
                                allkappa+=bands[j].heat_capacity[i]*bands[j].group_velocity[i]*bands[j].group_velocity[i]*bands[j].relaxation_time[i];
                                allheat+=bands[j].heat_capacity[i];
                            }
                            //cout<<endl;
                        }

                    }
                } else
                {
                    for (int j = 0; j < nband; ++j)
                    {
                        getline(inFile, line);
                        //cout<<line<<endl;
                        stringstream ss(line);

                        ss >> bands[j].group_velocity[i] >> bands[j].relaxation_time[i] >> bands[j].heat_capacity[i] >> bands[j].lattice_ratio[i] >> bands[j].heat_ratio[i];
                    }
                    break;
                }
                break;
                //cout<<endl;
            }
        }
    }

    /*getline(inFile, line);
    getline(inFile, line);
    bands.resize(Num_Matter);
    bands.resize(nband);
    for (int i = 0; i < nband; ++i)
    {
        bands[i].group_velocity.resize(Num_Matter);
        bands[i].relaxation_time.resize(Num_Matter);
        bands[i].heat_capacity.resize(Num_Matter);
        bands[i].lattice_ratio.resize(Num_Matter);
        bands[i].heat_ratio.resize(Num_Matter);
    }
    kappabulk.resize(Num_Matter);
    capacitybulk.resize(Num_Matter);
    for (int i = 0; i < Num_Matter; ++i)
    {
        while (getline(inFile, line))
        {
            //cout<<line.find("matter")<<" "<<line.length()<<endl;
            if (line.find("matter") > 0 && line.find("matter") < line.length())
            {
                //cout<<line<<endl;
                for (int j = 0; j < nband; ++j)
                {
                    getline(inFile, line);
                    //cout<<line<<endl;
                    stringstream ss(line);

                    ss >> bands[j].group_velocity[i] >> bands[j].relaxation_time[i] >> bands[j].heat_capacity[i] >> bands[j].lattice_ratio[i] >> bands[j].heat_ratio[i];
                }
                break;
            }

        }
    }*/
    for (int i = 0; i < Num_Matter; ++i)
    {
        kappabulk[i] = 0;
        capacitybulk[i]=0;
        for (int j = 0; j < nband; ++j)
        {
            if (Dimension_Material == 2)
                kappabulk[i] += bands[j].heat_capacity[i] * bands[j].group_velocity[i] * bands[j].group_velocity[i] * bands[j].relaxation_time[i] / 2;
            else if (Dimension_Material == 3)
                kappabulk[i] += bands[j].heat_capacity[i] * bands[j].group_velocity[i] * bands[j].group_velocity[i] * bands[j].relaxation_time[i] / 3;
            capacitybulk[i]+=bands[j].heat_capacity[i];
        }
        //cout<<kappabulk[i]<<endl;
        if(isnan(kappabulk[i]))
        {
            cout<<"Error: mistake in phonon properties, nan"<<endl;
            exit(0);
        }
        double sumL=0;
        double sumheat=0;
        for (int j = 0; j < nband; ++j) {
            sumL+= bands[j].lattice_ratio[i];
            sumheat+=bands[j].heat_ratio[i];
        }
        if(abs(sumL-1.0)>1e-6||abs(sumL-1.0)>1e-6)
        {
            cout<<"Error: mistake in phonon properties, enegy conservation does not obey"<<endl;
            exit(0);
        }
    }

    while (getline(inFile, line))
    {
        if (line.find("GeometryMaterialType") >= 0 && line.find("GeometryMaterialType") < line.length() && str.find(ant) == str.npos)
        {
            while (getline(inFile, line))
            {
                if (line.length() <= 1)
                {
                    break;
                }
                stringstream ss(line);
                int a, b;
                ss >> a >> b;
                if(abs(b)>100)
                {
                    break;
                }
                geo_matter_index.push_back(b-1);

            }
        }
    }
}
// BTEAngle::~BTEAngle()
// {
//     cout << "析构函数被调用" << endl;
// }
BTEBand::~BTEBand(){
    //cout << "BTEBand析构函数被调用" << endl;
}
