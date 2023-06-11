//
// Created by yuehu on 2023/5/21.
//

#include "TransientBTE/transient.h"
#include <algorithm>
#include <iomanip>
#include <chrono>
using namespace std;

void Transient::solve(int Use_Backup, double error_temp_limit, double error_flux_limit,double deltaT, double totalT,int use_TDTR,double pulse_time,double repetition_frequency,double modulation_frequency,double xy_r)
{
    double Num_Max_Iter=totalT/deltaT;
    _set_cell_matrix_larger();

    errorIncreaseTime=0;
    _set_initial(Use_Backup);
    
    
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++)
    {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local)
        {
            _get_bound_ee(iband_local,inf_local);
        }
    }
    _set_bound_ee_1();

    for (int ib = 0; ib < numBound; ++ib) {
        for (int inf_local = 0; inf_local < numDirectionLocal; ++inf_local) {
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
                int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
                for (int icell = 0; icell < 2; ++icell) {
                    int ie =boundaryCell[ib][icell];
                    int jface =boundaryFace[ib][icell];
                    if(ie>=0)
                    {
                        if (heatRatio[matter[ie]][iband][inf] != 0)
                        {
                            double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf]
                                                                                                                           * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                            if (dotproduct < 0)
                            {
                                if (boundaryType[ib]==1)
                                {
                                    double e = heatCapacity[matter[ie]][iband][inf] * boundaryThermal[ib];
                                    ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib*2] = e;
                                }
                                else if (boundaryType[ib]  == 2)
                                {
                                    double einsum1 = 0;
                                    double temp1 = 0;
                                    for (int nft = 0; nft < numDirectionLocal; ++nft)
                                    {
                                        double dotproduct1=(groupVelocityX[matter[ie]][iband][nft] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][nft]
                                                                                                                                        * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][nft] * elementFaceNormZ[jface + ie * 6]);
                                        if (dotproduct1 >= 0)
                                        {
                                            einsum1 += energyDensity[iband_local][nft][ib] *
                                                       (dotproduct1* modeWeight[matter[ie]][iband][nft]);
                                            temp1 += (dotproduct1*modeWeight[matter[ie]][iband][nft]);
                                        }
                                    }
                                    double e = einsum1 / temp1;
                                    ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib*2] = 0;//e;
                                }
                                else if (boundaryType[ib]  == 3)
                                {
                                    vec Reflectr;
                                    double dotproduct1=(directionX[inf] * elementFaceNormX[jface + ie * 6] + directionY[inf]
                                                                                                             * elementFaceNormY[jface + ie * 6] + directionZ[inf] * elementFaceNormZ[jface + ie * 6]);
                                    double ReflectrX=directionX[inf]- elementFaceNormX[jface + ie * 6] * dotproduct1 * 2;
                                    double ReflectrY=directionY[inf]- elementFaceNormY[jface + ie * 6] * dotproduct1 * 2;
                                    double ReflectrZ=directionZ[inf]- elementFaceNormZ[jface + ie * 6] * dotproduct1 * 2;
                                    //Reflectr = angles->direction[inf] - mesh->Elements[ie].faces[jface].norm * (angles->direction[inf] * mesh->Elements[ie].faces[jface].norm) * 2;
                                    double close = 1;
                                    int nf = -1;
                                    nf=0; //add ru
                                    for (int k = 0; k < numDirectionLocal; ++k)
                                    {
                                        double length= sqrt(pow(ReflectrX-directionX[k],2)+pow(ReflectrY-directionY[k],2)+pow(ReflectrZ-directionZ[k],2));
                                        if (length < close)
                                        {
                                            nf = k;
                                            close = length;
                                        }
                                    }
                                    ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib*2] = energyDensity[iband_local][nf][ib];

                                }

                                else if (boundaryType[ib] < 0)
                                {
                                    double e =energyDensity[iband_local][inf_local][boundaryCell[boundaryConnect[ib]][icell]];
                                    e = e + heatCapacity[matter[ie]][iband][inf] * (boundaryThermal[ib] - boundaryThermal[boundaryConnect[ib]]);

                                    ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib*2] = e;
                                }
                                else
                                {
                               }
                            }
                            else
                            {

                            }
                        }
                        else
                        {

                        }
                    }

                }
            }
        }
    }
    for (int iband = 0; iband < numBand; ++iband) {
        for (int inf = 0; inf < numDirection ; ++inf) {
            for (int ib = 0; ib < numBound*2; ++ib) {
               // cout<< ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib]<<endl;
            }

        }

    }

    auto total_iter_time = chrono::microseconds(0);
    auto get_gradient_time = chrono::microseconds(0);
    auto get_Re_time=chrono::microseconds(0);
    auto solver1_time = chrono::microseconds(0);
    auto Boundary_time = chrono::microseconds(0);
    auto set_vertex_time=chrono::microseconds(0);
    auto face_time=chrono::microseconds(0);
    auto non_frourier_time = chrono::microseconds(0);
    auto set_bound_time = chrono::microseconds(0);
    auto macro_bound_time = chrono::microseconds(0);
    auto macro_iter_time = chrono::microseconds(0);
    auto trasfer1_time = chrono::microseconds(0);

    _get_CellMatrix_larger();
    if(use_TDTR==1) {
        ofstream outputT("TTG.dat");
        ofstream outputTemptopave("TDTR.dat");
        for (int nt = 0; nt < Num_Max_Iter; ++nt) {
            total_iter_time = chrono::microseconds(0);
            get_gradient_time = chrono::microseconds(0);
            get_Re_time = chrono::microseconds(0);
            solver1_time = chrono::microseconds(0);
            Boundary_time = chrono::microseconds(0);
            set_vertex_time = chrono::microseconds(0);
            face_time = chrono::microseconds(0);
            non_frourier_time = chrono::microseconds(0);
            set_bound_time = chrono::microseconds(0);
            macro_iter_time = chrono::microseconds(0);
            trasfer1_time = chrono::microseconds(0);

            auto total_iter_start = chrono::high_resolution_clock::now();
            copy();
            for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
                for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
                    auto get_gradient_start = chrono::high_resolution_clock::now();
                    _get_gradient_larger(0, iband_local, inf_local);
                    auto get_gradient_end = chrono::high_resolution_clock::now();
                    get_gradient_time += chrono::duration_cast<chrono::microseconds>(
                            get_gradient_end - get_gradient_start);

                    auto get_Re_start = chrono::high_resolution_clock::now();
                    _get_explicit_Re(nt, 2, 0, iband_local, inf_local, deltaT);
                    auto get_Re_end = chrono::high_resolution_clock::now();
                    get_Re_time += chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                    auto solve_start = chrono::high_resolution_clock::now();

                    for (int icell = 0; icell < numCell; ++icell) {
                        //if(Re[icell]!=0)
                        //cout<<Re[icell]<<endl;
                        int inf = ((inf_local) * numProc + worldRank) % numDirection;
                        int iband =
                                iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
                        energyDensity[iband_local][inf_local][icell] = energyDensity[iband_local][inf_local][icell] *
                                                                       (1 - deltaT /
                                                                            relaxationTime[matter[icell]][iband][inf]) -
                                                                       deltaT * Re[icell];
                    }
                    auto solve_end = chrono::high_resolution_clock::now();
                    solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                    auto boundary_start = chrono::high_resolution_clock::now();
                    _get_bound_ee(iband_local, inf_local);
                    auto boundary_end = chrono::high_resolution_clock::now();
                    Boundary_time += chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                    _recover_temperature(iband_local, inf_local);
                    _get_total_energy(iband_local, inf_local);
                    _get_heat_flux(iband_local, inf_local);


                    MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                }
            }

            auto set_bound_start = chrono::high_resolution_clock::now();
            _set_bound_ee_1();
            auto set_bound_end = chrono::high_resolution_clock::now();
            set_bound_time += chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

            outputT << (nt + 1) * deltaT << " " << temperature[0] << endl;
            //memorize TDTR temp
            double Temptop_ave=0;
            double rr=0;
            double RR=xy_r*2;
            double savve=0;
            double heatratio=0;
            for(int icell=0;icell<numCell;icell++){
                rr = pow(elementCenterX[icell] - 0, 2) + pow(elementCenterY[icell ] - 0, 2);
                if(rr<pow(RR, 2))
                {
                    for(int jface=0;jface<elementFaceSize[icell];jface++){
                        if(elementFaceCenterZ[jface + icell * 6]==0) {
                            heatratio = exp(-2 * rr / pow(RR, 2));
                            savve = savve + heatratio;
                            Temptop_ave += temperature[icell] * heatratio;
                            break;
                        }
                    }
                }

            }
            Temptop_ave=Temptop_ave/savve;
            outputTemptopave << nt+1 << "   " << Temptop_ave <<endl;


            auto total_iter_end = chrono::high_resolution_clock::now();
            total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
            MPI_Barrier(MPI_COMM_WORLD);

            if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit)) {
                nt = Num_Max_Iter;
                if (worldRank == 0)
                    _print_out();
            }
            if (errorIncreaseTime >= 10000) {
                nt = Num_Max_Iter;
                if (worldRank == 0)
                    _print_out();
                if (worldRank == 0)
                    cout << "error increases for 10 times, maybe the solution can not converge. Try \"-1\" in limiter"
                         << endl;
                MPI_Barrier(MPI_COMM_WORLD);
            }
            if (nt % 1 == 0) {
                if (worldRank == 0)
                    _print_out();
            }

            if (worldRank == 0) {

                cout << "  Time taken by inner loop: " << 1.0 * total_iter_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by gradient 1: " << 1.0 * get_gradient_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by BTE Re: " << 1.0 * get_Re_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by BTE solver: " << 1.0 * solver1_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by Boundary: " << 1.0 * Boundary_time.count() / 1000 << " milliseconds"
                     << " " << 1.0 * set_bound_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by set_vertex: " << 1.0 * set_vertex_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by face " << 1.0 * face_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by non_fourier: " << 1.0 * non_frourier_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by macro: " << 1.0 * macro_iter_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by transfer: " << 1.0 * trasfer1_time.count() / 1000 << " milliseconds" << endl;
                cout << "----------------------------------------------------------------------------------" << endl;
            }

            if (worldRank == 0) {
                for (int j = 1; j < numProc; ++j) {
                    MPI_Send(&nt, 1, MPI_INT, j, 10, MPI_COMM_WORLD);
                }

            } else {
                MPI_Status status;
                MPI_Recv(&nt, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
            }
            //cout<<nt<<endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        outputT.close();
        outputTemptopave.close();
    }
    if(use_TDTR==0){
        ofstream outputT("TTG.dat");
        for (int nt = 0; nt < Num_Max_Iter; ++nt) {
            total_iter_time = chrono::microseconds(0);
            get_gradient_time = chrono::microseconds(0);
            get_Re_time = chrono::microseconds(0);
            solver1_time = chrono::microseconds(0);
            Boundary_time = chrono::microseconds(0);
            set_vertex_time = chrono::microseconds(0);
            face_time = chrono::microseconds(0);
            non_frourier_time = chrono::microseconds(0);
            set_bound_time = chrono::microseconds(0);
            macro_iter_time = chrono::microseconds(0);
            trasfer1_time = chrono::microseconds(0);

            auto total_iter_start = chrono::high_resolution_clock::now();
            copy();
            for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
                for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
                    auto get_gradient_start = chrono::high_resolution_clock::now();
                    _get_gradient_larger(0, iband_local, inf_local);
                    auto get_gradient_end = chrono::high_resolution_clock::now();
                    get_gradient_time += chrono::duration_cast<chrono::microseconds>(
                            get_gradient_end - get_gradient_start);

                    auto get_Re_start = chrono::high_resolution_clock::now();
                    _get_explicit_Re(nt, 2, 0, iband_local, inf_local, deltaT);
                    auto get_Re_end = chrono::high_resolution_clock::now();
                    get_Re_time += chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                    auto solve_start = chrono::high_resolution_clock::now();

                    for (int icell = 0; icell < numCell; ++icell) {
                        //if(Re[icell]!=0)
                        //cout<<Re[icell]<<endl;
                        int inf = ((inf_local) * numProc + worldRank) % numDirection;
                        int iband =
                                iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
                        energyDensity[iband_local][inf_local][icell] = energyDensity[iband_local][inf_local][icell] *
                                                                       (1 - deltaT /
                                                                            relaxationTime[matter[icell]][iband][inf]) -
                                                                       deltaT * Re[icell];
                    }
                    auto solve_end = chrono::high_resolution_clock::now();
                    solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                    auto boundary_start = chrono::high_resolution_clock::now();
                    _get_bound_ee(iband_local, inf_local);
                    auto boundary_end = chrono::high_resolution_clock::now();
                    Boundary_time += chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                    _recover_temperature(iband_local, inf_local);
                    _get_total_energy(iband_local, inf_local);
                    _get_heat_flux(iband_local, inf_local);


                    MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                }
            }

            auto set_bound_start = chrono::high_resolution_clock::now();
            _set_bound_ee_1();
            auto set_bound_end = chrono::high_resolution_clock::now();
            set_bound_time += chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

            outputT << (nt + 1) * deltaT << " " << temperature[0] << endl;
            auto total_iter_end = chrono::high_resolution_clock::now();
            total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
            MPI_Barrier(MPI_COMM_WORLD);

            if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit)) {
                nt = Num_Max_Iter;
                if (worldRank == 0)
                    _print_out();
            }
            if (errorIncreaseTime >= 10000) {
                nt = Num_Max_Iter;
                if (worldRank == 0)
                    _print_out();
                if (worldRank == 0)
                    cout << "error increases for 10 times, maybe the solution can not converge. Try \"-1\" in limiter"
                         << endl;
                MPI_Barrier(MPI_COMM_WORLD);
            }
            if (nt % 1 == 0) {
                if (worldRank == 0)
                    _print_out();
            }

            if (worldRank == 0) {

                cout << "  Time taken by inner loop: " << 1.0 * total_iter_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by gradient 1: " << 1.0 * get_gradient_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by BTE Re: " << 1.0 * get_Re_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by BTE solver: " << 1.0 * solver1_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by Boundary: " << 1.0 * Boundary_time.count() / 1000 << " milliseconds"
                     << " " << 1.0 * set_bound_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by set_vertex: " << 1.0 * set_vertex_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by face " << 1.0 * face_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by non_fourier: " << 1.0 * non_frourier_time.count() / 1000 << " milliseconds"
                     << endl;
                cout << "  Time taken by macro: " << 1.0 * macro_iter_time.count() / 1000 << " milliseconds" << endl;
                cout << "  Time taken by transfer: " << 1.0 * trasfer1_time.count() / 1000 << " milliseconds" << endl;
                cout << "----------------------------------------------------------------------------------" << endl;
            }

            if (worldRank == 0) {
                for (int j = 1; j < numProc; ++j) {
                    MPI_Send(&nt, 1, MPI_INT, j, 10, MPI_COMM_WORLD);
                }

            } else {
                MPI_Status status;
                MPI_Recv(&nt, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
            }
            //cout<<nt<<endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        outputT.close();
    }
        MPI_Barrier(MPI_COMM_WORLD);


    //_get_bound_temp();
    //_get_bound_flux();


#ifdef USE_TIME
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // cout << "end_iter1" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (worldRank == 0)
        _print_out();
    if (worldRank == 0)

#ifdef USE_TIME
        cout << "Time taken by iteration: " << duration.count() * 0.001 << " milliseconds" << endl;
#endif

    _delete_cell_matrix();


}

void Transient::solve_first_order(int Use_Backup, int Use_Limiter, double error_temp_limit,
                                  double error_flux_limit, double deltaT, double totalT) {

}


